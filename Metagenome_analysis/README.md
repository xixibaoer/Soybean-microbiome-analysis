# Soybean Rhizosphere Metagenome Analyses

**Author:** Qi Wang  
**Affiliation:** Max Planck Institute for Plant Breeding Research, Cologne, Germany

## Overview

This workflow analyzes soybean rhizosphere shotgun metagenomes to characterize microbial community composition and nutrient-associated functional genes across host genotypes and phosphate treatments. It combines Bash scripts for sequencing data processing with R code for statistical analysis and visualization. Computationally intensive steps include Slurm job configurations for an HPC cluster.

The primary design comprises 36 samples: two genotypes (WT and *phr1 phr2*), six phosphate treatments (0, 10, 20, 50, 100, and 200 µM), and three biological replicates per combination. Plotting sections also accommodate control samples where available.

## Workflow

| Section | Analysis | Main outputs |
| --- | --- | --- |
| 1 | Merge sequencing files within samples; fastp quality control; Bowtie2 host depletion | Filtered FASTQ files and QC reports |
| 2–3 | Kraken2 classification; Bracken abundance estimation and table merging | Taxonomic abundance matrices |
| 4–7 | Community ordination, composition, and diversity analysis | PCoA, PERMANOVA, family/genus plots, Shannon diversity |
| 8–9 | BBnorm coverage normalization; MEGAHIT co-assembly; QUAST assessment | Assembled contigs and assembly reports |
| 10–12 | Retrieve KO reference proteins and prepare sequence sets | KO-specific protein FASTA files |
| 13–14 | MAFFT alignment, trimAl trimming, HMM construction and indexing | Combined HMM database |
| 15–16 | Prodigal gene prediction, HMMER annotation, read mapping, abundance calculation, and KO assignment | Gene and KO-tagged abundance matrices |
| 17 | KO aggregation and comparison across treatments | Functional fold-change heatmaps |

Sections 11 and 12 provide alternative reference preparation approaches: CD-HIT clustering or retention of all downloaded reference sequences.

## Requirements

### Software

- **Read processing:** fastp, Bowtie2, SAMtools.
- **Taxonomic profiling:** Kraken2, Bracken, and its `combine_bracken_outputs.py` utility.
- **Assembly:** BBTools (`bbnorm.sh`), MEGAHIT, QUAST.
- **Functional annotation:** curl, CD-HIT, MAFFT, trimAl, Prodigal, HMMER, SeqKit.
- **Execution:** Bash, standard Unix utilities, Python 3 with pandas, and Slurm for HPC submission.
- **R analysis:** vegan, ggplot2, readr, dplyr, tibble, readxl, tidyverse, tidyr, scales, and patchwork.

Rendering the analysis document as HTML additionally requires R Markdown tooling, Pandoc, and rmdformats.

### Input data and references

1. Paired-end FASTQ files, organized by biological sample.
2. An Excel metadata file with a `Sheet1` worksheet containing `SampleID`, `Genotype`, and `Treatment`.
3. A soybean Bowtie2 reference index.
4. A Kraken2 database with matching Bracken read-length distribution files.
5. A text file listing selected KEGG Orthology identifiers, one per line.
6. A functional metadata table containing `KO`, `Gene`, and `Nutrient` for heatmap annotation.

Example biological sample identifiers are `WT_0uM_1` and `phr1_2_0uM_1`. Raw-read discovery uses `*_R1_*.fq.gz` and `*_R2_*.fq.gz`. Downstream scripts use paired filenames such as:

```text
WT_0uM_1_R1_nohost.fq.gz
WT_0uM_1_R2_nohost.fq.gz
```

Match sample identifiers to the metadata and use consistent filenames between processing stages. Treatment labels in the plotting code include `Ctl` and `Brady_*`.

## Running the analysis

1. Replace `/Data/to/path` placeholders with the appropriate input, output, and reference locations.
2. Update Conda activation commands, Slurm settings, and log directories for the computing environment.
3. Run read preprocessing in Section 1, followed by taxonomic profiling and community analysis in Sections 2–7.
4. Run coverage normalization and co-assembly in Sections 8–9. Prepare the functional reference database using Sections 10–14.
5. Complete gene annotation and abundance estimation in Sections 15–16, then generate functional heatmaps with Section 17.

Bash sections can be extracted into shell scripts and submitted to Slurm; R sections can be run in R or extracted into R scripts. For example, using your chosen filenames:

```bash
sbatch analysis_step.sh
Rscript community_analysis.R
```

These filenames are illustrative. The source contains both executable and `eval=FALSE` chunks; select the sections needed for the analysis. Resource settings are specified per job and should reflect the available cluster resources.

## Analysis settings

| Step | Settings used in the script |
| --- | --- |
| Taxonomic profiling | Kraken2 confidence 0.1; Bracken read length 150 bp; species through phylum levels |
| Community ordination | Mean taxon abundance >0.01%; Hellinger transformation; Bray–Curtis distance; PCoA |
| Statistical analysis | PERMANOVA with 999 permutations; dispersion assessment; Shannon diversity comparisons |
| Coverage normalization | BBnorm `target=100`, `min=5` |
| Assembly | MEGAHIT minimum contig length 1,000 bp |
| Functional references | Up to 200 proteins per KO; optional CD-HIT clustering at 99% identity |
| HMM profiles | MAFFT, trimAl, and hmmbuild; minimum five reference sequences |
| Functional search | HMMER sequence reporting threshold `-E 1e-10` |

Functional gene abundance uses length-normalized mapped counts and within-sample TPM scaling:

```text
RPK = mapped count / gene length in kilobases
TPM = RPK / sum of RPK values in the sample × 1,000,000
```

TPM is relative to the gene set included in the calculation. These DNA-based measurements represent functional gene abundance rather than transcriptional expression. Section 17 aggregates abundances by KO and calculates fold changes relative to the WT 0 µM reference.

## Outputs and reproducibility

The workflow produces quality-control reports, host-depleted reads, taxonomic abundance matrices, community plots and statistics, assembled contigs, reference HMMs, predicted genes, annotation tables, mapping/count files, and functional heatmaps.

Retain sample metadata, software versions, database releases, reference genome information, execution parameters, and job logs with the final results. Sequencing data accession numbers, data availability details, and a repository license can be added when established.
