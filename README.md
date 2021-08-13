---

<h3>Introduction</h3>  

---  

Here we present all necessary code to reproduce the results (incl. Figures and Tables) of Key et al. 2021. We use [`miniconda`](https://conda.io/en/latest/miniconda.html) for package management, [`snakemake`](https://snakemake.readthedocs.io/en/stable/) for pipeline processing, and `python3` for analysis and visualisation, with few bits of `matlab` and `R` left within. We are utilizing a ton of python modules all included in the spyder4 conda environment. Thank you py-community! While I aim to make the analyses as accessible as possible, a medium level of bioinformatic expertise is required for succesful execution.

All snakemake pipelines have a similar structure, made of the `Snakefile`, `cluster.slurm.json` (config requires adjustment for individual cluster environment), `snakemakeslurm.sh` (requires minor individual adjustments for execution of snakemake on HPC), `logs` for `stdout`/`stderr` output, and optional an `envs/` and `scripts/` folder containing conda environment files to set up the various bioinformatic tools utilized or custom-made scripts, respectively. All snakemakes are designed for a high-performance-computing cluster using the submission system `slurm`. `Snakefile`'s need to be revised for individual usage. Each snakemake analysis requires in addition to all provided files/folders a `samples.csv` input file. An example `samples.csv` can be found within each respective snakemake folder. The `samples.csv` contains information about individual single isolates: 
- `Path` describes the folder path including the bgzip raw sequencing data (incl. trailing `/`)
- `Sample` is the isolate name
- `ReferenceGenome` the reference genome identifier
- `ProviderName` the unique raw data file name (excluding `_R1.fastq.gz`/`_R2.fastq.gz`)
- `Subject` patient identifier

The raw genomic sequence data for 1,336 S. aureus isolates has been deposited at SRA Bioproject [PRJNA715375](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715375/) and [PRJNA715649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715649/). Details about each isolate, incl. individual Biosample SRR-IDs, are presented in Supplementary Table 10 of the publication, which we also provide as csv here `SOMTable10.tsv`. The entire raw sequence data should be obtained and stored prior to the analysis (eg.  using `fastq-dump`).

The pipeline is divided into the following major analyses:
1. Within-patient
2. Across-patient
3. MGE
4. Generation of Figures and Tables

Below details are provided about the execution of each step.

We hope you find the code useful. In case you utilize it for your own analyses, please cite: ` On-person adaptive evolution of Staphylococcus aureus during atopic dermatitis increases disease severity
Felix M. Key, Veda D. Khadka, Carolina Romo-González, Kimbria J. Blake, Liwen Deng, Tucker C. Lynn, Jean C. Lee, Isaac M. Chiu, Maria Teresa García-Romero, Tami D. Lieberman
bioRxiv 2021.03.24.436824; doi: https://doi.org/10.1101/2021.03.24.436824 `

---

<h3>1. Within-patient analysis</h3>  

---

The analysis splits into two parts. First, we use four snakemake scripts to generate the patient-specific pangenomes, perform basic data QC, alignment, and variant calling. Second, within-patient analyses incl. variant and isolate filtering, and phylogenetic, adaptive evolution and molecular clock analysis.

<h3>Raw data processing</h3>  

Snakemake processing splits into four parts:
1. Basic filtering and taxonomic classification (`snakemake/withinpat/kraken2`) 
 - The kraken2 database has been made with all refseq genomes (archaea bacteria plasmid viral human fungi protozoa UniVec) following default recommendations.
 - Execution (similar for all snakemake runs):
`bash run_snakemake.slurm`
2. Pangenome assembly and annotation (`snakemake/withinpat/assembly`)
 - Builds upon identified S.aureus isolates (step 1) for generating pangenome
 - Python3 script utilize the following modules: `argparse`,`sys`,`subprocess`,`gzip`,`os`,`SeqIO` from `Bio`,`statistics`
 - [`spades.py`](https://github.com/ablab/spades) executable required, which has to be specified in the `Snakefile`
 - 
4. Alignment and variant calling (`snakemake/withinpat/mapping`)
5. Building candidate_mutation_table for each patient (`snakemake/withinpat/case`)





To be continued.



