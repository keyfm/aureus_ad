---

<h3>Introduction</h3>  

---  

Here we present all necessary to reproduce the results (incl. Figures and Tables) of Key et al. 2021. We use [`snakemake`](https://snakemake.readthedocs.io/en/stable/) for raw data processing and `python3` for analysis and visualisation, with few bits of `matlab` and `R` left. The raw genomic sequence data for 1,336 S. aureus isolates has been deposited at SRA Bioproject [PRJNA715375](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715375/) and [PRJNA715649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715649/). Details about each isolate, incl. individual Biosample SRR-IDs, are presented in Supplementary Table 10 of the publication, which we also provide as csv here `SOMTable10.tsv`. The entire raw sequence data should be obtained and stored prior to the analysis (eg.  using `fastq-dump`).

The pipeline is divided into four major analyses:
1. Within-patient
2. Across-patient
3. MGE analysis
4. Generation of Figures and Tables

Below details are provided about the execution of each step.

We hope you find the code useful. In case you utilize major parts of it for your own analyses, we would appreciate if you acknowledge our publication. `XXX`

---

<h3>Within-patient analysis</h3>  

---

The analysis splits into two parts. First, we use four snakemake scripts to generate the patient-specific pangenomes, perform basic data QC, alignment, and variant calling. Second, within-patient analyses incl. variant and isolate filtering, and phylogenetic, adaptive evolution and molecular clock analysis.

<h3>Within-patient - data processing</h3>  

Snakemake processing splits into four parts:
1. Basic filtering and species assignment
2. Pangenome assembly and annotation
3. Variant calling
4. Building candidate_mutation_table

To be continued.


Come back soon!

