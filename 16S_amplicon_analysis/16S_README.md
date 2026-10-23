# Soybean Rhizosphere Microbiome: 16S Amplicon Analysis

**Author:** Qi Wang  
**Affiliation:** Max Planck Institute for Plant Breeding Research, Cologne, Germany

## Overview

This workflow analyzes soybean-associated bacterial communities using 16S rRNA gene amplicon sequence variant (ASV) tables. It compares community structure, alpha diversity, and taxonomic composition across soil types, sample compartments, and microbial inoculation treatments.

The analysis combines R Markdown, phyloseq, vegan, and ggplot2. It starts from processed ASV counts, taxonomy assignments, sample metadata, and representative sequences. Raw FASTQ processing and ASV inference are performed upstream and are not included in this script.

## Experimental factors

| Factor | Categories used in the analysis |
| --- | --- |
| Soil type | Clay, Loamy, Sandy |
| Sample compartment | Soil, Rhizosphere, Endosphere |
| Inoculation treatment | Control, Bradyrhizobium, SC3, SC4, SC6 |

The alpha-diversity section standardizes `SC3`, `SC4`, and `SC6` to `SynCom3`, `SynCom4`, and `SynCom6`. Sample metadata determine the actual sample numbers and replication.

## Workflow

| Section | Analysis | Main outputs |
| --- | --- | --- |
| 1–2 | Load packages and define input paths | Configured analysis environment |
| 3 | Import ASV counts, taxonomy, metadata, and representative sequences | Separate phyloseq objects for each compartment |
| 4–5 | Merge compartment datasets, remove empty entries, and inspect object structure | Combined phyloseq object and RDS file |
| 6 | Calculate Bray–Curtis distances, perform PCoA, and run PERMANOVA | Soil-specific ordination plots and test results |
| 7 | Calculate Observed ASV richness and Shannon diversity; run ANOVA and Tukey comparisons separately within each soil | Diversity tables, statistical models, and annotated plots |
| 8 | Visualize family-level taxonomic composition | Relative abundance plots |

## Requirements

Use an R installation compatible with the required CRAN and Bioconductor packages.

**Main analysis packages:** phyloseq, Biostrings, vegan, ggplot2, tidyverse, dplyr, tibble, viridis, multcompView, and patchwork.

**Additional packages loaded by the supplied script:** dada2, ggrepel, DESeq2, ggforce, gridExtra, data.table, ggpubr, ggtree, plotly, ComplexHeatmap, scales, ggsci, reshape2, and colorspace.

Rendering the HTML report requires rmarkdown, knitr, rmdformats, and Pandoc. The family-composition section also requires the external `graphical_methods.R` file defining `plot_composition()`.

Record the R version, package versions, and taxonomy database release used for the final analysis.

## Input files

Provide the following files for each compartment, using the `RS`, `Root`, and `Soil` entries in `data_paths`:

| Input | Expected structure |
| --- | --- |
| ASV count table | Tab-separated file; `OTU.ID` column followed by sample count columns |
| Taxonomy table | Tab-separated file; `Feature ID` and a semicolon-separated `Taxon` column |
| Sample metadata | Tab-separated file with `SampleID`, `Soil_Type`, and `Treatment` |
| Representative sequences | FASTA file with ASV identifiers matching the count and taxonomy tables |

The script assigns compartment labels using the dataset of origin: `RS` becomes `Rhizosphere`, `Root` becomes `Endosphere`, and `Soil` becomes `Soil`.

Sample identifiers must be unique across compartments and correspond to the metadata. ASV identifiers must represent the same sequences wherever they are shared between datasets. Diversity calculations require nonnegative count data.

## Running the analysis

1. Replace `/Data/to/path` placeholders with the input, output, and reference paths for your environment.
2. Provide the processed ASV tables, taxonomy, metadata, and representative sequences for the three compartments.
3. Run the import and phyloseq construction sections before the community analyses.
4. Run the soil-specific ordination sections and the alpha-diversity section.
5. Configure the composition-plot inputs and the external plotting helper before running Section 8.

Run sections interactively in RStudio or render the configured `.Rmd` document with `rmarkdown::render()`. Review chunk execution and caching options for the intended run.

## Community structure

Section 6 calculates Bray–Curtis distances and uses PCoA to visualize community structure separately for Clay, Loamy, and Sandy soils. Points are colored by inoculation treatment and shaped by compartment.

The supplied models test compartment (`Origin`) in Clay and Loamy soils and treatment (`Treatment`) in Sandy soil, using 9,999 permutations. Each plot identifies the factor tested and reports its PERMANOVA p-value and R².

## Alpha diversity

Section 7 calculates **Observed ASV richness** and **Shannon diversity** from `meta_seq_clean` and joins the measurements to metadata by `SampleID`.

The default compartment is:

```r
target_origin <- "Endosphere"
```

Change this to `"Rhizosphere"` or `"Soil"` to analyze another compartment.

For each diversity metric, the script fits a separate one-way ANOVA within each soil type:

```r
Observed ~ Treatment
Shannon ~ Treatment
```

Tukey-adjusted pairwise comparisons are summarized using compact letter displays. Letters apply only within the corresponding soil and metric: treatments sharing a letter are not significantly different at the 0.05 threshold under that test. Letters are not comparisons between soils.

The models assume independent observations and appropriate residual behavior. Observed richness should be interpreted alongside sequencing depth and sampling coverage.

## Taxonomic composition

Section 8 uses the external `plot_composition()` function to display family-level relative abundances, with the top 20 families and plots faceted by soil type. Its helper file and plotting inputs must be provided with the analysis.

## Outputs

- Combined phyloseq dataset: `merged_and_filtered_phyloseq.rds`.
- Soil-specific PCoA plots with PERMANOVA annotations.
- Alpha-diversity measurements: `Alpha_diversity_<compartment>.csv`.
- Per-soil ANOVA models, Tukey comparisons, and letter annotations: `Alpha_diversity_statistics_<compartment>.rds`.
- Combined richness and Shannon plots: `Alpha_diversity_<compartment>.pdf`.
- Family-level community composition plots.

Output locations are defined in the relevant sections. Retain input metadata, upstream ASV-processing methods, reference database information, and `sessionInfo()` with the final results.
