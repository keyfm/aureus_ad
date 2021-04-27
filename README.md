---

<h3>Introduction</h3>  

---  

Here we present all necessary to reproduce the results (incl. Figures and Tables) of Key et al. 2021. We use [`miniconda`](https://conda.io/en/latest/miniconda.html) for package management, [`snakemake`](https://snakemake.readthedocs.io/en/stable/) for pipeline processing, and `python3` for analysis and visualisation, with few bits of `matlab` and `R` left. We are utilizing a ton of python modules all included in the spyder4 conda environment. Thank you py-community! While I aim to make the analyses as accessible as possible, a medium level of bioinformatic expertise is required for succesful execution.

All snakemake pipelines have a common structure, made of the `Snakefile`, `cluster.slurm.json` (config requires adjustment for individual cluster environment), `snakemakeslurm.sh` (requires minor individual adjustments for execution of snakemake on HPC), `logs` for `stdout`/`stderr` output, and optional an `envs/` and `scripts/` folder containing conda environment files to set up the various bioinformatic tools utilized or custom-made scripts, repectively. All snakemakes are designed for a high-performance-computing cluster using the submission system `slurm`. Each snakemake analysis requires in addition to all provided files/folders a `samples.csv` input file. An example `samples.csv` can be found within each respective snakemake folder. The `samples.csv` contains all/some of the following items: 
- `Path` describes the folder path including the bgzip raw sequencing data (incl. trailing `/`)
- `Sample` is the isolate name
- `ReferenceGenome` the reference genome identifier
- `ProviderName` the unique raw data file name (excluding `_1.fastq.gz`/`_2.fastq.gz`)
- `Subject` patient identifier

The raw genomic sequence data for 1,336 S. aureus isolates has been deposited at SRA Bioproject [PRJNA715375](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715375/) and [PRJNA715649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715649/). Details about each isolate, incl. individual Biosample SRR-IDs, are presented in Supplementary Table 10 of the publication, which we also provide as csv here `SOMTable10.tsv`. The entire raw sequence data should be obtained and stored prior to the analysis (eg.  using `fastq-dump`).

The pipeline is divided into the following major analyses:
1. Within-patient
2. Across-patient
3. MGE
4. Generation of Figures and Tables

Below details are provided about the execution of each step.

We hope you find the code useful. In case you utilize it for your own analyses, please cite: `XXX`

---

<h3>Within-patient analysis</h3>  

---

The analysis splits into two parts. First, we use four snakemake scripts to generate the patient-specific pangenomes, perform basic data QC, alignment, and variant calling. Second, within-patient analyses incl. variant and isolate filtering, and phylogenetic, adaptive evolution and molecular clock analysis.

<h3>Raw data processing</h3>  

Snakemake processing splits into four parts:
1. Basic filtering and taxonomic classification (`snakemake/withinpat/kraken2`) 
 - The kraken2 database has been made with all refseq genomes (archaea bacteria plasmid viral human fungi protozoa UniVec) following default recommendations.
 - Execution:
`sbatch --mem=2000 -o logs/0master.log -c 1 --job-name='SM.master' --time 2-00:00:00 -p defq,sched_mem1TB_centos7 --wrap="\
 bash snakemakeslurm.sh ; "`
2. Pangenome assembly and annotation (`snakemake/withinpat/assembly`)
3. Variant calling (`snakemake/withinpat/align`)
4. Building candidate_mutation_table for each patient (`snakemake/withinpat/case`)

 sbatch --mem=2000 -o logs/0master.log -c 1 --job-name='SM.master' --time 2-00:00:00 -p defq,sched_mem1TB_centos7 --wrap="\
 bash snakemakeslurm.sh ; "





To be continued.



