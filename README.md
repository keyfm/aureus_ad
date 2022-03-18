
<h3>S. aureus evolution on atopic dermatitis skin</h3>  

---

<h3>Introduction</h3>  

---  

Here we present all necessary code to reproduce the results (incl. major Figures and Tables) of Key et al. 2021. We use [`miniconda`](https://conda.io/en/latest/miniconda.html) for package management, [`snakemake`](https://snakemake.readthedocs.io/en/stable/) for pipeline processing, and `python3` for analysis and visualisation, with few bits of `matlab` and `R` left within. We are utilizing a ton of python modules all included in the provided spyder4 conda environment. Thank you py-community! While I aim to make the analyses as accessible as possible, a medium level of bioinformatic expertise is required for succesful execution.

All snakemake pipelines have the same structure, made of the `Snakefile`, `cluster.slurm.json` (config requires adjustment for individual cluster environment), `snakemakeslurm.sh` (requires minor individual adjustments for execution of snakemake on HPC), `run_snakemake.slurm`, `logs` for stdout/stderr output, and optional an `envs` and `scripts` folder containing conda environment files to set up the various bioinformatic tools utilized and custom-made scripts, respectively. All snakemakes are designed for a high-performance-computing cluster using the submission system `slurm`. `Snakefile`'s need to be revised for individual usage (ie. paths). Each snakemake analysis requires in addition to all provided files/folders a `samples.csv` input file. An example `samples.csv` can be found within each respective snakemake folder. The `samples.csv` contains information about individual single isolates: 
- `Path` describes the folder path including the bgzip raw sequencing data (incl. trailing `/`)
- `Sample` is the isolate name
- `ReferenceGenome` the reference genome identifier (folder label that contains `genome.fasta(.gz)` file)
- `ProviderName` the unique raw data file name (excluding `_R1.fastq.gz`/`_R2.fastq.gz`)
- `Subject` patient identifier

The raw genomic sequence data for 1,336 S. aureus isolates has been deposited at SRA Bioproject [PRJNA715375](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715375/) and [PRJNA715649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715649/). Details about each isolate, incl. individual Biosample SRR-IDs, are presented in Supplementary Table 10 of the publication, which we also provide as csv here `SOMTable10.tsv`. The entire raw sequence data should be obtained and stored prior to the analysis (eg.  using `fastq-dump`).

The pipeline is divided into the following major analyses:
1. Across-patient
2. Within-patient
4. Public data analysis
5. Generation of Figures and Tables

Below details are provided about the execution of each step.

We hope you find the code useful. In case you utilize it for your own analyses, please cite: `On-person adaptive evolution of Staphylococcus aureus during atopic dermatitis increases disease severity`
Felix M. Key, Veda D. Khadka, Carolina Romo-González, Kimbria J. Blake, Liwen Deng, Tucker C. Lynn, Jean C. Lee, Isaac M. Chiu, Maria Teresa García-Romero, Tami D. Lieberman
bioRxiv 2021.03.24.436824; doi: https://doi.org/10.1101/2021.03.24.436824


---

<h3>1. Across-patient analysis</h3>  

---

The analysis splits into two parts. First, we use two snakemake scripts to perform basic data QC, alignment to a reference genome, and variant calling.  Second, the across-patient analyses incl. variant and isolate filtering, and the generation of all input for RAxML to build the phylogeny. The entire analysis is designed to run on a HPC.

<h3>Raw data processing</h3>  

Snakemake processing splits into two parts:
1. `mapping`: Alignment and variant calling using the genomic data of all isolates.
2. `case`: Building one `candidate_mutation_table.pickle.gz` for all patients. 

<h3>Analysis</h3>  

 Filter variants and build a multi-fasta file as input for a tool of your choice to  build a maximum-likelihood phylogeny. Set up the `spyder4_full_env.yml` environment. Run `across_patient_analysis.py` within that conda environment. 


---

<h3>2. Within-patient analysis</h3>  

---

The analysis splits into two parts. First, we use four snakemake scripts to generate the lineage-specific pangenomes, perform basic data QC, alignment, and variant calling. They are designed to run on an HPC. Second, within-patient analyses incl. variant and isolate filtering, and phylogenetic, adaptive evolution and molecular clock analysis. That analysis is designed to run on a regular laptop/desktop.

<h3>Raw data processing</h3>  

Snakemake processing splits into four parts:
1. `kraken2`: Basic filtering and taxonomic classification
 - The kraken2 database has been made with all refseq genomes (archaea bacteria plasmid viral human fungi protozoa UniVec) following default recommendations.
2. `assembly`: Pangenome assembly and annotation 
 - Builds upon identified S.aureus isolates (step 1) for generating pangenome
 - Python3 script utilize the following modules: `argparse`,`sys`,`subprocess`,`gzip`,`os`,`SeqIO` from `Bio`,`statistics`
 - [`spades.py`](https://github.com/ablab/spades) executable required, which has to be specified in the `Snakefile`
4. `mapping`: Alignment and variant calling using genomic data of all isolates.
5. `case`: Building `candidate_mutation_table.pickle.gz` for each patient. Has to be run for each patient (ie. assembled genome) individually within the case folder (`case/subject_XX/`). 

<h3>Analysis</h3>  

 Here we produce all results of the within-person evolution analysis, utilizing `candidate_mutation_table.pickle.gz` for each patient. Set up the `spyder4_full_env.yml` environment. Run `within_patient_analysis.py` within the provided conda environment `spyder4_full_env.yml`. 


---

<h3>3. Public Data Analysis</h3>  

---

In order to reproduce the public data analysis first download from SRA the raw read info for each isolate shown in Extended Data Table 9. For each isolate, build the and annotate an assembly using the snakemake `assembly` provided in the `within_host_analysis`. Next, blast the query sequence fasta (Extended Data Table 8) against each assembly (with the following flags `-outfmt 5 -max_hsps 1`). The assembled genomes and the annotation (`gff` format) as well as the `xml` blast output are the input for the `public_data_analysis.py` (see within script).

---

<h3>4. Figures and Tables </h3>  

---

The `figures_and_tables_generator.py` contains all the code necessary to rebuild the figures and tables presented in the publication. Additional metadata is provided in the `metadata` folder if necessary. 

The evolvograms (aka muller plots) are generated using a mix of custom-python code, as well as py/R packages. The python-wrapper `xxx` and the required conda environment `zzz` is available.


If there are any questions or comments please reach out to fkey [at] mit . edu


