"""Repair the invalid normal-math run pattern emitted by Pandoc OMML."""

from pathlib import Path
from zipfile import ZipFile

INVALID_NORMAL_RUN = b'<m:rPr><m:nor/><m:sty m:val="p"/></m:rPr>'
VALID_NORMAL_RUN = b"<m:rPr><m:nor/></m:rPr>"


def repair_pandoc_omml(docx_path: Path) -> int:
    """Remove Pandoc's redundant style element from normal math runs."""
    docx_path = Path(docx_path)
    with ZipFile(docx_path) as source:
        entries = [(info, source.read(info.filename)) for info in source.infolist()]

    replacements = 0
    repaired_entries = []
    for info, data in entries:
        if info.filename == "word/document.xml":
            replacements = data.count(INVALID_NORMAL_RUN)
            data = data.replace(INVALID_NORMAL_RUN, VALID_NORMAL_RUN)
        repaired_entries.append((info, data))

    if replacements == 0:
        return 0

    temporary_path = docx_path.with_suffix(".repairing.docx")
    with ZipFile(temporary_path, "w") as target:
        for info, data in repaired_entries:
            target.writestr(info, data)
    temporary_path.replace(docx_path)
    return replacements
