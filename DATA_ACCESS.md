# Data access and release boundary

## GitHub materials

The public GitHub repository may contain analysis code, manuscript files,
documentation, figures, tables, and aggregate statistical results. It must not
contain participant-level recordings or linkage information.

## Borealis publication dataset

The Human Research Ethics Board of the University of Victoria gave ethical
approval for this work under protocol H22-00451. The study team confirmed
authorization to publish the de-identified participant-level data prepared for
this manuscript. Borealis has reserved DOI
[`10.5683/SP4/7YSHYF`](https://doi.org/10.5683/SP4/7YSHYF).

The prepared Borealis package contains the data required to reproduce the
reported heel-IMU analysis:

- 32 participant archives under codes `ID_NNN`
- 128 bilateral heel IMU CSV files
- 126 synchronization pulse CSV files
- 897 cleaned GaitRite workbooks
- A data dictionary, inventory, checksums, and reproduction documentation

## Excluded materials

Do not upload these materials to GitHub or Borealis:

- Sacrum or hip IMU recordings
- Source participant identifiers, initials, or collection dates
- Source folder names or code-linkage files
- TDMS or HDF5 source files
- Unredacted GaitRite workbooks
- Files outside the validated publication package

The publication copy retains only the two heel-sensor blocks, rebases absolute
timestamps, removes identifying workbook content, and uses `ID_NNN` codes. The
validated package reproduces all 10,472 matched stride pairs used for the
frozen analysis. Never modify the source recordings.
