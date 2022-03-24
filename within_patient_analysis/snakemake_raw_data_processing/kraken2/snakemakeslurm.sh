#!/bin/bash

# launch snakemake to run jobs via SLURM
SM_PARAMS="job-name ntasks partition time mail-user mail-type error output"

SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb}"

for P in ${SM_PARAMS}; do SM_ARGS="$SM_ARGS --$P {cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"

# our SLURM error/output paths expect a logs/ subdir in PWD
mkdir -p logs


		
### run snakemake
snakemake -p \
    $* \
     --latency-wait 10 \
    -j 1000 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch $SM_ARGS" \
    --cluster-status scripts/slurm_status.py \
    --restart-times 2\
	--keep-going \
	--nolock \
	--reason \
	--use-conda \
	--max-jobs-per-second 100 \
	--max-status-checks-per-second 1000\
	--restart-times 3 \
	--rerun-incomplete \
	--conda-prefix /scratch/mit_lieberman/tools/conda_snakemake/

    # --dag \
    # | dot -Tsvg > dag.svg
