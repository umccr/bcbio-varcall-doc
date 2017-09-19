#!/bin/bash
#SBATCH -n 1
#SBATCH -p vccc
#SBATCH --mem=8000
#SBATCH --time=24:00:00
#SBATCH -o log.nortum_varscan.out
#SBATCH -e log.nortum_varscan.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pdiakumis@unimelb.edu.au

bcbio_nextgen.py \
    ../config/cancer-normal.yaml \
    -n 64 \
    -q vccc \
    -s slurm \
    -t ipython \
    -r 'timelimit=1-00:00:00' \
    --retries 1 \
    --timeout 120
