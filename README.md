<h3>Introduction</h3>  

---  

Here we present all necessary to reproduce the results (incl. Figures and Tables) of Key et al. 2021. For the analysis we use primarily [`snakemake`](https://snakemake.readthedocs.io/en/stable/) and `python3`, with few bits of `matlab` and `R` left. The raw genomic sequence data for 1336 S. aureus isolates has been deposited at SRA Bioproject [PRJNA715375](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715375/) and [PRJNA715649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715649/). Details about each isolate, incl. individual Biosample SRR-IDs, are presented in Supplementary Table 10 of the publication, which we also provide as csv here `SOMTable10.tsv`. The entire raw data sequence data should be obtained and stored prior to the analysis (eg.  using `fastq-dump`).

The pipeline is divided into four major analyses:
1. Within-patient
2. Across-patient
3. MGE analysis
4. Generation of Figures and Tables

Below details are provided about the execution of each step.

We hope you find the code presented useful. In case you utilize major parts of it for your own analyses, we would appreciate if you cite our publication.

<h3>Within-patient analysis</h3>  

---

The analysis utilizes patient-specific pangenomes for alignment 

Site is currently under construction.

Come back soon!

