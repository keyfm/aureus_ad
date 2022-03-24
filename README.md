
# _S. aureus_ evolution on atopic dermatitis skin </h1>  

Here we present the computational code part of the analysis presented in: 

Felix M. Key, Veda D. Khadka, Carolina Romo-González, Kimbria J. Blake, Liwen Deng, Tucker C. Lynn, Jean C. Lee, Isaac M. Chiu, Maria Teresa García-Romero, Tami D. Lieberman
bioRxiv 2021.03.24.436824; doi: https://doi.org/10.1101/2021.03.24.436824

We hope you find the code useful. In case you recycle it for your own analyses please cite our study.



<h2>Introduction</h2>  


The analysis uses [`miniconda`](https://conda.io/en/latest/miniconda.html) for package management, [`snakemake`](https://snakemake.readthedocs.io/en/stable/) for pipeline processing, and `python3` for analysis and visualisation, with few bits of `matlab` and `R` left within. We are utilizing multiple python modules all included in the provided conda environment (`spyder4_full_env.yml`). Thank you py-community! While I aim to make the analyses as accessible as possible, a medium level of bioinformatic expertise is required for succesful execution.

The analyses is divided into four parts:
1. Across-patient analysis
2. Within-patient analysis
3. Metaanalysis Public Data
4. Figure and Table generation

Each part (except the latter) starts with raw data processing done via snakemake, followed by the analysis in python3.

All snakemake pipelines have the same structure, made of the `Snakefile`, `cluster.slurm.json` (config file requires adjustment for individual cluster environment), `snakemakeslurm.sh` (requires minor individual adjustments for execution of snakemake on HPC), `run_snakemake.slurm`, `logs` (for stdout/stderr output), and (sometimes) an `envs` and `scripts` folder containing conda environment files and custom-made scripts, respectively. All snakemakes are designed for execution on a computing cluster using the submission system `slurm`. `Snakefile`'s need to be revised for individual usage (ie. paths). Information about each isolate dataset is fed into Snakemake via a `samples.csv` input file. An example file is provided in each snakemake folder. The `samples.csv` contains the following information about individual single isolates: 
- `Path` is the path with the raw sequencing data stored (incl. trailing `/`)
- `Sample` is the isolate name
- `ReferenceGenome` the reference genome identifier (folder label that contains `genome.fasta(.gz)` file)
- `ProviderName` the unique raw data file name within `Path` (excluding `_R1.fastq.gz`/`_R2.fastq.gz`)
- `Subject` patient identifier

The raw genomic sequence data for 1,587 *S. aureus* isolates is available from SRA Bioproject [PRJNA715375](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715375/), [PRJNA715649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA715649/) and [PRJNA816913](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA816913/). Details about each isolate, incl. individual Biosample SRR-IDs, are presented in the publication (Supplementary Table 10). The entire raw sequence data should be obtained and stored prior to the analysis (eg.  using `fastq-dump`).



<h2>1. Across-patient analysis</h2>  


The analysis splits into two parts. First, we use two snakemake scripts to perform basic data QC, alignment to a reference genome, and variant calling.  Second, the across-patient analyses incl. variant and isolate filtering, and the generation of all input for RAxML to build the phylogeny. The entire analysis is designed to run on a HPC.

<h4>Raw data processing</h4>  

Source: `across_patient_analysis/snakemake_raw_data_processing`
1. `mapping`: Alignment to reference genome and variant calling using the genomic data of all isolates.
2. `case`: Building one `candidate_mutation_table.pickle.gz` for all patients that contains basecalls and relevant summary statisitcs for downstream analysis. 

<h4>Analysis</h4>  

Source: `across_patient_analysis/across_patient_analysis.py`

Using `candidate_mutation_table.pickle.gz` to filter variants and build a multi-fasta file as input phylogeny reconstruction. Set up the `spyder4_full_env.yml` environment. Run `across_patient_analysis.py` within that conda environment. 



<h2>2. Within-patient analysis</h2>  


Again, the analysis splits into two parts. First, we use four snakemake scripts to generate the lineage-specific pangenomes, perform basic data QC, alignment, and variant calling. They are designed to run on an HPC. Second, within-patient analyses incl. variant and isolate filtering, and phylogenetic, adaptive evolution and molecular clock analysis. That analysis is designed to run on a regular laptop/desktop.

<h4>Raw data processing</h4>

Source: `within_patient_analysis/snakemake_raw_data_processing`

1. `kraken2`: Basic filtering and taxonomic classification
 - The kraken2 database has been made with all refseq genomes (archaea bacteria plasmid viral human fungi protozoa UniVec) following default recommendations.
2. `assembly`: Pangenome assembly and annotation 
 - Builds upon identified *S.aureus* isolates (step 1) for generating pangenome
 - Python3 script utilize the following modules: `argparse`,`sys`,`subprocess`,`gzip`,`os`,`SeqIO` from `Bio`,`statistics`
 - [`spades.py`](https://github.com/ablab/spades) executable required, which has to be specified in the `Snakefile`
4. `mapping`: Alignment and variant calling using genomic data of all isolates.
5. `case`: Building `candidate_mutation_table.pickle.gz` for each patient. Has to be run for each patient (ie. assembled genome) individually within the case folder (`case/subject_XX/`). 

<h4>Analysis</h4>  

Source: `within_patient_analysis/within_patient_analysis.py`

 Here we produce all results of the within-person evolution analysis, utilizing individual `candidate_mutation_table.pickle.gz` for each patient. Run `within_patient_analysis.py` within the provided conda environment `spyder4_full_env.yml`. 



<h2>3. Metaanalysis Public Data</h2>  

Source: `metaanalysis_public_data/public_data_analysis.py`

In order to reproduce the public data analysis first download from SRA the raw read info for each isolate shown in Extended Data Table 9. For each isolate, build and annotate an assembly using the snakemake `assembly` provided in the `within_host_analysis`. Next, blast the query sequence fasta (Extended Data Table 8) against each assembly (with the following flags `-outfmt 5 -max_hsps 1`). The assembled genomes and the annotation (`gff` format) as well as the `xml` blast output are the input for the `public_data_analysis.py`.



<h2>4. Figures and Tables </h2>  

Source: `figures_tables/figures_and_tables_generator.py`

The `figures_and_tables_generator.py` contains all the code necessary to rebuild the figures and tables presented in the publication. Required metadata is provided in the `metadata` folder. 

<h4>Muller plots aka evolvograms </h4>  

Source: `muller_evolvograms/muller_plots.py`

The muller plots (also called evolvograms) are generated using the the tool [`lolipop`](https://github.com/cdeitrick/Lolipop) and the R-package [`ggmuller`](https://github.com/robjohnnoble/ggmuller), which are wrapped together within the custom python code in `muller_plots.py`. 




