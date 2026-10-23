# Soybean Microbiome Analysis

Scripts and documentation for analyzing soybean-associated microbial communities using **16S rRNA gene amplicon sequencing** and **shotgun metagenomics**. The workflows examine microbial community structure, diversity, and functional potential in relation to soil conditions, microbial inoculation, and host phosphate-starvation signaling.

## Repository contents

| Folder | Description |
| --- | --- |
| [16S_amplicon_analysis](16S_amplicon_analysis/) | ASV-table processing, community ordination, alpha diversity, soil-specific ANOVA/Tukey comparisons, and taxonomic composition plots. |
| [Metagenome_analysis](Metagenome_analysis/) | Read quality control, host depletion, taxonomic profiling, assembly, functional annotation, and gene abundance analysis. |

Each folder contains analysis scripts and its own README describing inputs, dependencies, and execution.

## Analyses

**16S amplicon analysis** compares Soil, Rhizosphere, and Endosphere communities across Clay, Loamy, and Sandy soils and microbial inoculation treatments. It starts from processed ASV counts, taxonomy, representative sequences, and sample metadata.

**Metagenome analysis** examines community composition and nutrient-associated functional genes across WT and *phr1 phr2* soybean plants exposed to different phosphate treatments. It includes Kraken2/Bracken profiling, MEGAHIT assembly, and HMM-based functional annotation.

## Getting started

1. Open the relevant analysis folder and read its README.
2. Install the required R packages and command-line software.
3. Prepare the input data and references, and update file paths and computing settings.
4. Run the analysis sections in the documented order.

Raw sequencing data and reference databases are supplied separately. Software requirements differ between workflows; computationally intensive metagenome steps include Slurm configurations for HPC execution.

## Author

**Qi Wang**  
Max Planck Institute for Plant Breeding Research, Cologne, Germany  
Contact: qwang2@mpipz.mpg.de
