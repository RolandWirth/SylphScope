<img width="1164" height="368" alt="image" src="https://github.com/user-attachments/assets/c2c522bf-684d-4926-a1ee-99dd11107907" />

# SylphScope - Sylph MPA Explorer

Interactive R Shiny application for exploratory analysis, visualisation and statistical comparison of Sylph MPA taxonomic abundance tables.

## Main features

- Upload one or multiple Sylph MPA tables (`.sylphmpa`, `.txt`, `.tsv`)
- Optional sample metadata with `SampleID` and `Group`
- Sample-depth normalisation of `sequence_abundance`
- MPA-style species/strain filtering to reduce double counting
- Phylum composition plots
- Species heatmaps
- Bray-Curtis PCoA and Hellinger PCA
- Bray-Curtis dendrogram
- PERMANOVA
- LEfSe-like biomarker screening
- ANI = 100 species/strain summaries
- Export of processed tables, plots, HTML report, citation metadata and session information

## Input format

The app expects Sylph MPA-style tables containing at least:

- `clade_name`
- `sequence_abundance`

Optional columns:

- `ANI`
- `Coverage`

Metadata file must contain:

- `SampleID`
- `Group`

## Normalisation

For each sample, the app calculates a clean sample depth after MPA-style strain filtering. Sequence abundance values are then normalised as:

```text
normalised_reads = sequence_abundance × (average_sample_sum / sample_sum)
normalised_percent = (normalised_reads / average_sample_sum) × 100
```

## Run locally

```r
shiny::runApp()
```

or from a terminal:

```bash
Rscript app.R
```

## Installation notes

The app contains an automatic package checker/installer. If you prefer manual installation, set this line in `app.R`:

```r
AUTO_INSTALL_PACKAGES <- FALSE
```

Then install the required packages manually from CRAN/Bioconductor.

## Citation

If you use this software, please cite:

> Wirth R. SylphScope: an interactive Shiny application for exploratory analysis of Sylph MPA taxonomic profiles. Manuscript in preparation.

## License

MIT License.
