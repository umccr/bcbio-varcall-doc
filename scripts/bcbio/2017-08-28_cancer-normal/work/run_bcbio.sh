#!/bin/bash
#SBATCH -n 1
#SBATCH -p vccc
#SBATCH --mem=8000
#SBATCH --time=240:00:00
#SBATCH -o log.nortum1.bcbio.out
#SBATCH -e log.nortum1.bcbio.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pdiakumis@unimelb.edu.au

bcbio_nextgen.py \
    ../config/cancer-normal.yaml \
    -n 64 \
    -q vccc \
    -s slurm \
    -t ipython \
    -r 'timelimit=5-00:00:00' \
    --retries 1 \
    --timeout 120
