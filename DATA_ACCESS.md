# Data access and release boundary

## Public materials

The public GitHub repository may contain analysis code, manuscript files,
documentation, figures, tables, and aggregate statistical results. The
aggregate Borealis deposit uses reserved DOI `10.5683/SP4/7YSHYF`.

## Restricted materials

The University of Victoria Human Research Ethics Board approved the study
under protocol H22-00451. The approved protocol restricts access to
de-identified data to approved research personnel and states that data will
not be kept on the web.

Do not upload any of the following to GitHub, Borealis, or another web service
without written UVic REB clarification or an approved amendment:

- Raw or processed heel IMU recordings
- Sacrum or hip IMU recordings
- GaitRite trial exports
- GR Pulse or synchronization records
- Participant-coded stride, walk, or event-level outputs
- Collection dates, initials, source folder names, or linkage files

## Offline preparation

A future participant-level release would require new ethics authorization.
After authorization, create a separate publication copy. Remove the sacrum
sensor block, retain and relabel the two heel sensors, replace source names
with approved study codes, rebase absolute timestamps, and validate the copy
against the frozen 10,472-pair analysis. Never modify the source recordings.
