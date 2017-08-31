Snakemake
=========

### Activate snakemake environment
```
source activate snake
```

### Run snakemake on Spartan
```
snakemake -j 2 --cluster-config cluster.json \
  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time}" \
  -s snake.snakefile
```
