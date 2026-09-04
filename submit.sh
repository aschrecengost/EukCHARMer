#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH -p cpu,uri-cpu,cpu-preempt
#SBATCH --mem=64g
# module load uri/main
# module load conda/latest
# module load snakemake/6.10.0-foss-2021b 
snakemake \
-s Snakefile --jobs 100 --use-conda --conda-frontend conda \
--cluster "sbatch -p {cluster.partition} --job-name={rule} --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.ntasks} --nodes={cluster.nodes} --output=slurm-logs/%x.%j.log --verbose" \
--cluster-config cluster.yaml --latency-wait 120