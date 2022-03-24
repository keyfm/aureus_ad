#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=1-00:00:00
#SBATCH -o 00log_stderrout.txt
#SBATCH -e 00log_stderrout.txt
#SBATCH --mem=8000
#SBATCH --mail-user=USERNAME@DOMAIN
#SBATCH --mail-type=FAIL
snakemake --unlock
bash snakemakeslurm.sh 
