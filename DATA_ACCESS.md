# Data access

The public dataset is hosted by **Borealis, UVic Research Data Collection**.

- DOI [10.5683/SP4/7YSHYF](https://doi.org/10.5683/SP4/7YSHYF)
- Version 1.0, released September 17, 2026
- 128 heel IMU CSVs, 126 synchronization CSVs, 897 GaitRite XLSX workbooks
  and five supporting files

The deposit uses 32 anonymized study codes, `ID_NNN`, in the three data folders
shown in [the format guide](docs/data-format.md). Inventory archive fields
refer to the participant ZIPs used to prepare the deposit, not the layout
required by the pipeline.

`ID_186` and `ID_195` lack L2 synchronization files, so those conditions cannot
enter paired analysis. `ID_191` was excluded because bilateral heel recordings
were unavailable. The deposited exclusions table documents these omissions.
No other participant-condition record is excluded based on error or pair count.

Verify downloads using the deposited `SHA256SUMS_EXTRACTED_DATA.txt`. Its hashes
cover 1,151 raw files, not supporting CSV tables converted to TAB by Borealis.

## Data boundary

This archive contains analysis code only.
The small screening-group table in `metadata/` transcribes already reported
Supplementary Table S1 values and preserves its bootstrap ordering.
Participant recordings and generated participant-level outputs belong outside
Git and are not release assets.

Published data retain the two heel-sensor blocks, use relative timestamps and
omit source identifiers, linkage information and collection dates. They do not
include sacrum recordings, HDF5/TDMS originals or unredacted workbooks. Do not
add these materials to this repository.

Study ethics approval is H22-00451. Data reuse follows the dataset's repository
terms. The code's MIT licence does not relicense the data.
