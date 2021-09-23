#!/bin/bash
#SBATCH -q batch -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-5:00
#SBATCH --mem=16GB
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR
module load singularity
export G=https://raw.githubusercontent.com/TheJacksonLaboratory/mousegwas
nextflow run TheJacksonLaboratory/mousegwas -profile jax --yaml $G/example/gait_nowild_noBL.yaml --shufyaml $G/example/gait_shuffle_noBL.yaml --input $G/example/gait_paper_strain_survey_2019_11_12.csv --outdir gait_output_noBL --clusters 1 --addpostp "--colorgroup --meanvariance --set3 --minherit 0.25 --loddrop 1.5" --addheatmap "--meanvariance -p 0.1"  -resume
