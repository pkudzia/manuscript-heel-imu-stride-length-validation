# Paper 1 final citation audit

Date: 2026-08-26

## Scope and status

- Audited every citation-supported claim in the final manuscript against the cited source.
- Checked all 35 references rendered in the revised manuscript.
- Also checked Nilsson et al. (2012), which the prior draft cited for an unsupported orientation-error claim.
- Used the PDFs in `literature/` where available. For sources without a local PDF, used the original article text or authoritative publication record.
- Made only claim-level corrections needed to match the sources. No numerical study result changed.
- Promoted the matched audited Word, LaTeX, and PDF set after author approval. The prior originals are preserved in `archive/2026-08-26_pre_audit_and_superseded/backups/`.

## Recommended citation decisions

| Decision | Recommendation | Rationale |
|---|---|---|
| Add Webster et al. (2005) and Arens et al. (2021) to the criterion-system sentence | Accept | Webster validates GaitRite against motion capture. Arens uses optical motion capture as the criterion. Both already appear in the bibliography. |
| Move Foxlin (2005) to the fixed-threshold and tuning claim | Accept | Foxlin supports detector tuning and stationary-state assumptions. It does not support the draft's Parkinson-specific example. |
| Replace Nilsson et al. (2012) and Mahony et al. (2008) with Foxlin (2005) and Hannink et al. (2017) for the orientation-error claim | Accept | Nilsson describes an open-source ZUPT-aided INS. Mahony describes the orientation filter. Foxlin explains how tilt error creates horizontal acceleration error, and Hannink documents the foot-trajectory consequence and swing-phase accelerometer gating. |
| Remove Lin et al. (2025) from the stops-and-turns sentence and retain Mariani et al. (2010) only for prior turn testing | Accept | Mariani included turning tasks. Lin is a review of sensor placement and does not establish a stop-or-turn advantage. |

No new paper is required. These changes only reassign references already in the bibliography. Nilsson et al. (2012) becomes uncited and therefore does not render in the revised reference list.

## Changes made in the versioned audited set

| Location | Prior wording or implication | Audited wording or action | Reason |
|---|---|---|---|
| Introduction, opening | Combined mobility decline, falls, and cognitive deterioration under one broad topic sentence | Split the evidence into variability predicting falls, slower gait predicting falls and cognitive decline, and stride length at or below 0.64 m predicting adverse events and mortality | The cited studies support different outcomes. The revision prevents each citation from appearing to support every outcome. |
| Introduction, criterion systems | Instrumented walkways and motion capture measure stride length accurately, without citations | Added Webster et al. (2005) and Arens et al. (2021) | Both sources directly document the criterion systems. |
| Introduction, sensor factors | Assigned hardware noise, attachment, placement, and algorithm claims jointly to Wahlstrom and Kuderle | Assigned accumulated sensor error to Wahlstrom and placement or attachment effects to Kuderle | Each source now supports the specific claim attached to it. |
| Introduction, threshold detector | Claimed the detector fails during Parkinsonian shuffling and cited Foxlin | Removed the Parkinson-specific example. Stated that fixed settings can fail when movement differs from the tuning data | Foxlin supports tuning sensitivity, not that clinical example. |
| Introduction, bilateral ZUPT | Presented the paper's full acceptance window and recovery procedure as if Arens described it | Limited the Arens-supported claim to contralateral toe-off identification. Marked the acceptance window and recovery search as this study's implementation | Arens supports the bilateral biomechanical relationship, not every parameter used here. |
| Introduction, Sijobert study | Characterised prior validation as steady walking in healthy adults | Specified 10 healthy adults and 12 people with Parkinson's disease during straight, single-task walks | This matches the reported cohort and protocol. |
| Introduction, dual task | Said dual-task walking slows gait and increases variability | Changed to “can slow gait and increase stride-to-stride variability” | The cited studies show task- and cohort-dependent effects, not a universal response. |
| Introduction, Boutaayamou study | Described one dual-task condition but omitted the exact cohort and outcome scope | Specified intra-session reliability in 101 middle-aged adults and no criterion validation | This matches the study design. |
| Introduction, Laidig study | Said pathological cohorts walked at fixed speeds | Changed to self-selected speed on an instrumented treadmill | The original study used self-selected speeds. |
| Introduction, clinical relevance | Claimed dual-task assessment supports cognitive monitoring | Narrowed this to fall-risk evaluation and gait-cognition interactions | The cited work supports association and screening relevance, not established monitoring performance. |
| Introduction, heel evidence | Called Kuderle the only large heel-specific study and described variability | Specified a 14-person multi-position study with higher heel single-stride error, especially during fast walking | This accurately describes its sample, comparison, and error outcome. |
| Introduction, Arens study | Described wearable-mounted sensors and bilateral validation without the reference details | Specified bilateral foot sensors, contralateral toe-off, eight people post-stroke, and optical motion capture | These details match the source. |
| Methods, stance Step 1 | Attached the cited methods to the study-specific 40 ms window | Cited the signal features and stated separately that this implementation used 40 ms | Neither cited source establishes this study's exact window. |
| Methods, stance Step 2 | Attached Arens to the exact 150 ms acceptance window | Cited Arens for contralateral toe-off confirmation and stated separately that this implementation used 150 ms | Arens does not specify this implementation's exact window. |
| Methods, gravity residual | Cited Nilsson and Mahony for horizontal residual acceleration caused by orientation error | Cited Foxlin and Hannink and changed “inflate” to “bias” | The new citations directly support the mechanism and consequence. Bias is direction-neutral. |
| Methods, agreement statistics | Applied Bland and Altman plus Shrout and Fleiss to a mixed list of metrics | Assigned bias and limits of agreement to Bland and Altman, and ICC(2,1) to Shrout and Fleiss | The revision matches the scope of each statistical source. |
| Discussion, Kuderle comparison | Called the heel result greater single-stride variability | Changed it to higher single-stride error | Kuderle reports estimation error across placements, not biological stride variability. |
| Discussion, Ensink comparison | Said Ensink found no accuracy decrement under irregular stepping | Specified sensors on both feet, sternum, and lower back, healthy and stroke cohorts, and only small changes in method disagreement | This avoids overstating a non-significant or absolute null effect. |
| Discussion, turns | Suggested the bilateral detector's intended turn advantage and cited Mariani plus Lin | Stated that prior foot-worn systems included turns, then stated that this detector remains untested during stops and turns | Neither source validates an advantage of this bilateral detector. |
| Discussion, sensor locations | Described Trojaniello and Washabaugh as foot- and trunk-mounted studies | Changed the locations to feet, ankles, and shanks | Those are the sensor locations used by the cited studies. |

## Source-by-source audit matrix

| Source | Claim checked | Audit result |
|---|---|---|
| Arens et al. (2021) | Bilateral foot sensors, contralateral toe-off, post-stroke cohort of eight, optical motion-capture reference | Supported after separating the paper's 150 ms window and recovery logic from the cited method. |
| Beauchet et al. (2005) | Cognitive dual-task effects on gait speed and variability | Supported with “can” rather than a universal causal statement. |
| Bland and Altman (1986) | Bias and 95% limits of agreement | Supported after separating these measures from MAE, ICC, and proportional bias. |
| Boutaayamou et al. (2025) | Heel IMU, 101 middle-aged adults, cognitive dual task, intra-session reliability, no criterion validation | Supported after specifying the study scope. |
| Bytyci et al. (2021) | Baseline stride length at or below 0.64 m predicts adverse events and mortality | Supported. |
| Commandeur et al. (2018) | Dual-task gait measures and fall-risk assessment in older adults | Supported after narrowing “cognitive monitoring” to fall-risk and gait-cognition relevance. |
| Ensink et al. (2023) | Healthy and stroke cohorts, irregular stepping, sensors on both feet, sternum, and lower back | Supported after replacing “no decrement” with the more precise finding of small changes in method disagreement. |
| Foxlin (2005) | Foot-mounted ZUPT, stationary updates, threshold tuning, and tilt-error effects on horizontal acceleration | Supported. Removed the unsupported Parkinson-specific example. |
| Greene et al. (2010) | Signal-derived adaptation of gyroscope-based gait-event thresholds | Supported. The manuscript identifies its own numerical thresholds as study-specific. |
| Hannink et al. (2017) | Foot-trajectory estimation, orientation-error consequences, and suppressing accelerometer orientation updates during swing | Supported. The manuscript describes a similar movement-gating rationale, not an identical gain schedule. |
| Hausdorff et al. (2008) | Cognitive dual-task effects on gait variability | Supported with non-universal wording. |
| Kottner et al. (2011) | GRRAS guidance for reporting reliability and agreement | Supported. |
| Küderle et al. (2022) | Fourteen healthy young adults, six foot locations, higher heel single-stride error, strongest difference during fast walking | Supported after correcting “variability” to “error.” |
| Laidig et al. (2021) | Calibration-free foot IMU validation in pathological cohorts at self-selected treadmill speed | Supported after correcting “fixed” to “self-selected.” |
| Li et al. (2023) | Adaptive ZUPT thresholds across walking-speed changes | Supported. |
| Lin et al. (2025) | Heel placement is uncommon in the reviewed IMU gait literature | Supported. Removed from the unsupported stops-and-turns claim. |
| Luo et al. (2024) | Heel-mounted validation showed greater spread than an embedded in-shoe sensor while mean errors remained small | Supported. |
| Mahony et al. (2008) | Nonlinear complementary attitude filter using gyroscope and accelerometer information | Supported for the filter description. Removed from the separate orientation-residual claim. |
| Maki (1997) | Stride-to-stride variability predicts falls | Supported. |
| Mariani et al. (2010) | Foot-worn IMU stride-length validation and inclusion of turning tasks | Supported. The manuscript no longer implies that it validates this bilateral detector. |
| Montero-Odasso et al. (2012) | Slow gait predicts cognitive decline and dual-task gait has gait-cognition relevance | Supported. |
| Morris et al. (2019) | APDM Opal spatiotemporal gait validity in older or clinical populations | Supported. |
| Nasreddine et al. (2005) | Montreal Cognitive Assessment basis and conventional score interpretation | Supported. The manuscript clearly labels its groups as study-defined screening groups rather than diagnoses. |
| Rampp et al. (2015) | Low-centimetre stride-length error from a lower-leg or foot-adjacent IMU in geriatric patients | Supported. |
| Salarian et al. (2004) | Sagittal-plane gyroscope basis for gait-event detection | Supported. The manuscript states that it adapted the method. |
| Shrout and Fleiss (1979) | ICC(2,1) framework | Supported after assigning this citation only to ICC. |
| Sijobert et al. (2015) | Straight, single-task validation in 10 healthy adults and 12 people with Parkinson's disease using foot and shank sensors | Supported after correcting the cohort description. |
| Skog et al. (2010) | Fixed-threshold zero-velocity detection using gyroscope magnitude and accelerometer variability | Supported. |
| Springer et al. (2006) | Dual-task walking and fall-risk relevance in older adults | Supported under the narrowed clinical-relevance wording. |
| Trojaniello et al. (2014) | Shank-mounted IMU gait and stride validation in older and pathological cohorts | Supported. Corrected the location summary in the Discussion. |
| Verghese et al. (2009) | Slow gait and gait variability as fall predictors | Supported after separating these outcomes from cognitive decline. |
| Wahlström et al. (2021) | Sensor error accumulation, ZUPT principles, and sensitivity to detector settings and walking speed | Supported. The BibTeX key retains the early-access year `Wahlstrom2020`, while the rendered year is 2021. |
| Washabaugh et al. (2017) | APDM Opal validity and repeatability with sensors on feet and ankles | Supported. Corrected the location summary in the Discussion. |
| Webster et al. (2005) | GaitRite concurrent validity against motion capture, with ICC 0.92 to 0.99 across speed, cadence, step length, and step time | Supported. The manuscript does not imply that the full range applies only to step length. |
| Yogev-Seligmann et al. (2008) | Cognitive dual-task effects on gait | Supported with non-universal wording. |
| Nilsson et al. (2012), removed from rendered references | Open-source foot-mounted ZUPT-aided INS | Does not support the former orientation-residual claim. No other final claim requires this citation. |

## Verification

- Both canonical LaTeX copies have identical SHA-1 hashes.
- The audited `manuscript_submission.docx` preserves the audited text and passes package validation against the user's edited Word file.
- The audited LaTeX compiles through BibTeX and three PDFLaTeX passes with no undefined citations or references.
- The audited PDF contains 19 letter-sized pages.
- Visual inspection of all 19 pages found no clipped figures, split tables, or interrupted sentences.
- Nilsson et al. (2012) does not render in the audited reference list.
