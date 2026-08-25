#!/usr/bin/env python3
"""Produce citations/verified.jsonl audit trail from references.bib.

For each entry with a DOI, hit CrossRef and lock in title, authors, year,
journal, and final URL. Writes one JSONL record per entry with
verification_status, source, and today's date.
"""

import json
import re
import sys
import time
import unicodedata
import urllib.request
from datetime import date
from pathlib import Path

import bibtexparser


def normalize(s: str) -> str:
    """Lowercase, strip accents, collapse all unicode hyphens to ascii hyphen."""
    s = s.lower()
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    for h in "‐‑‒–—―−":
        s = s.replace(h, "-")
    return s


BIB = Path(__file__).parent.parent / "references.bib"
OUT = Path(__file__).parent / "verified.jsonl"
TODAY = date.today().isoformat()


def latex_clean(s: str) -> str:
    s = re.sub(r"\\['`\"~^]?\{?([a-zA-Z])\}?", r"\1", s)
    s = re.sub(r"\\[a-zA-Z]+\s*\{?([^}]*)\}?", r"\1", s)
    s = s.replace("{", "").replace("}", "").replace("--", "-")
    return re.sub(r"\s+", " ", s).strip()


def crossref_lookup(doi: str, timeout=15):
    url = f"https://api.crossref.org/works/{doi}"
    try:
        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": (
                    "IMU_GaitSync-audit/1.0 (mailto:pkudzia@trackyourcore.com)"
                )
            },
        )
        with urllib.request.urlopen(req, timeout=timeout) as r:
            if r.status != 200:
                return None, f"HTTP {r.status}"
            data = json.loads(r.read().decode("utf-8"))
        return data.get("message"), None
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code}"
    except Exception as e:
        return None, str(e)


def main():
    with BIB.open() as fh:
        db = bibtexparser.load(fh)
    print(f"Parsed {len(db.entries)} entries from {BIB.name}", file=sys.stderr)

    records = []
    for entry in db.entries:
        key = entry.get("ID", "?")
        etype = entry.get("ENTRYTYPE", "?")
        bib_title = latex_clean(entry.get("title", ""))
        bib_authors = latex_clean(entry.get("author", entry.get("editor", "")))
        bib_year = entry.get("year", "").strip()
        bib_journal = latex_clean(entry.get("journal", entry.get("booktitle", "")))
        bib_doi = entry.get("doi", "").strip()
        bib_url = f"https://doi.org/{bib_doi}" if bib_doi else ""

        rec = {
            "bibkey": key,
            "entry_type": etype,
            "title": bib_title,
            "authors": bib_authors,
            "year": bib_year,
            "journal": bib_journal,
            "doi": bib_doi,
            "url": bib_url,
            "date": TODAY,
        }

        if not bib_doi:
            rec["verification_status"] = "needs_manual_verification"
            rec["notes"] = "no DOI in bib (book or other); verify via WorldCat/ISBN"
            records.append(rec)
            print(f"  [{key}] no DOI - manual", file=sys.stderr)
            continue

        msg, err = crossref_lookup(bib_doi)
        time.sleep(0.1)
        if err or not msg:
            rec["verification_status"] = "verify_failed"
            rec["error"] = err or "no message"
            rec["notes"] = "DOI does not resolve in CrossRef; entry may be fabricated"
            records.append(rec)
            print(f"  [{key}] FAIL: {err}", file=sys.stderr)
            continue

        cr_title = (msg.get("title") or [""])[0]
        cr_authors_list = msg.get("author") or []
        cr_journal = (msg.get("container-title") or [""])[0]
        cr_year_dp = msg.get("published-print", msg.get("published-online", {})).get(
            "date-parts"
        )
        cr_year = str(cr_year_dp[0][0]) if cr_year_dp else ""
        cr_volume = str(msg.get("volume", ""))
        cr_pages = msg.get("page", msg.get("article-number", ""))

        rec["crossref_title"] = cr_title
        rec["crossref_authors_first3"] = [
            f"{a.get('family', '?')}, {a.get('given', '?')}"
            for a in cr_authors_list[:3]
        ]
        rec["crossref_journal"] = cr_journal
        rec["crossref_year"] = cr_year
        rec["crossref_volume"] = cr_volume
        rec["crossref_pages"] = cr_pages

        bib_t = re.sub(r"\W", "", normalize(bib_title))
        cr_t = re.sub(r"\W", "", normalize(cr_title))
        bib_t_short = bib_t[:40]
        cr_t_short = cr_t[:40]
        title_match = bool(bib_t) and (bib_t_short in cr_t or cr_t_short in bib_t)

        bib_first = (
            normalize(bib_authors.split(" and ")[0].split(",")[0].strip())
            if bib_authors
            else ""
        )
        cr_first = (
            normalize(cr_authors_list[0].get("family", "")) if cr_authors_list else ""
        )
        author_match = bool(bib_first) and (
            bib_first == cr_first or bib_first in cr_first or cr_first in bib_first
        )

        rec["title_match"] = title_match
        rec["author_match"] = author_match

        if title_match and author_match:
            rec["verification_status"] = "verified_online"
            rec["notes"] = f"verified against CrossRef on {TODAY}"
        else:
            rec["verification_status"] = "verify_failed"
            rec["notes"] = (
                f"title_match={title_match}, author_match={author_match}; "
                f"bib first='{bib_first}' crossref='{cr_first}'"
            )

        records.append(rec)
        flag = "OK" if rec["verification_status"] == "verified_online" else "MISMATCH"
        print(f"  [{key}] {flag}", file=sys.stderr)

    with OUT.open("w") as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    n_ok = sum(1 for r in records if r["verification_status"] == "verified_online")
    n_manual = sum(
        1 for r in records if r["verification_status"] == "needs_manual_verification"
    )
    n_fail = sum(1 for r in records if r["verification_status"] == "verify_failed")
    print(
        (
            f"\n{len(records)} records: {n_ok} verified, "
            f"{n_manual} manual, {n_fail} failed"
        ),
        file=sys.stderr,
    )
    print(f"Wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
