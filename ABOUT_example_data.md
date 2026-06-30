# SylphScope corrected complex synthetic example data

This dataset contains 12 synthetic `.sylphmpa` files and matching metadata.

## Contents

- `example_data/`: 12 `.sylphmpa` files
- `metadata.csv`: metadata for Shiny upload
- `metadata.tsv`: same metadata as TSV
- `dataset_summary.tsv`: check table showing each sample sums to 100%

## Test features

- Four groups: Control, Treatment_A, Treatment_B, Treatment_C
- Three replicates per group
- Uneven sequencing depths
- Bacteria and Archaea
- Multiple phyla
- Species-level and strain-level terminal taxa
- Several ANI = 100 entries
- Group-specific synthetic abundance shifts for PCA, heatmaps, and biomarker-style testing

Upload all `.sylphmpa` files together, then upload `metadata.csv` or with `metadata.tsv`.
