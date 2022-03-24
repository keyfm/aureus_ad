#!/bin/bash

#
# check for clean
#
# https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-do-i-remove-all-files-created-by-snakemake-i-e-like-make-clean
# if [ "$1" == "clean" ]; then
#     echo 'rm $(snakemake --summary | tail -n+2 | cut -f1)'
#     snakemake --summary | tail -n+2 | cut -f1
#     rm -f $(snakemake --summary | tail -n+2 | cut -f1)
#     exit 0
# fi

#
# launch snakemake to run jobs via SLURM
# ntasks
SM_PARAMS="job-name ntasks partition time mail-user mail-type error output"

SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb}"

for P in ${SM_PARAMS}; do SM_ARGS="$SM_ARGS --$P {cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"

# our SLURM error/output paths expect a logs/ subdir in PWD
mkdir -p logs

### run snakemake
# -j defines total number of jobs executed in parallel on cluster
# -n dryrun
# -p print command lines
snakemake -p \
    $* \
    --use-conda \
    --latency-wait 30 \
    -j 600 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch $SM_ARGS" \
    --cluster-status /scratch/mit_lieberman/scripts//slurm_status.py \
    --rerun-incomplete \
    --restart-times 3 \
    --keep-going \
    --conda-prefix /scratch/mit_lieberman/tools/conda_snakemake/


    # --dag \
    # | dot -Tsvg > dag.svg
