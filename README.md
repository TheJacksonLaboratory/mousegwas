MouseGWAS
=========

Introduction
------------

This package was built to manage the GWAS analysis of mouse phenotypes. The mice in the study were genotypes using either MDA or UCLA chips and deposited in the mouse phenome database (https://phenome.jax.org/genotypes). 

Installation
------------
```bash
Rscript -e 'library(devtools); install_github("TheJacksonLaboratory/mousegwas")'
```

Input
-----

The input for the script is the genotype csv files downloaded from the MPD website, the measured phenotypes as a csv file and a yaml file describing the input file.
The input csv file should contain a column for strain, a column for sex and columns for phenotype measurements. The names of the columns should be defined in the yaml file using the keywords `strain` and `sex` and the phenotypes should be a list under the `phenotypes` keyword.
Another data that should reside in the yaml file is translation of strains to the strain names in the genotypes files, it is a dictionary under the `translate` keyword, and `F1` keyword which is a dictionary translating the F1 names to their parent names, make sure the female parent is always first, it will be used to determine the X chromosome of make F1s. Confounding SNPs could be given using the `confSNPs`, this might be useful to control for obvious markers like coat color alleles. For sanity check you can supply coat color under `coat` as a dictionary from strain name to coat color and execute a GWAS of coat color with `--coat_phenotype`, it can also be used as a covariate with `--coat_covar`.

Execution
---------

The script `run_GWAS.py` will read the input files and will prepapre the input for either GEMMA or pyLMM. In the case of GEMMA it will download version 0.98 to the working directory if it can't find the GEMMA executable, if you wish to use pyLMM it should be installed and available in the path. A common process would be creating a virtual environment in python, activating it and installing pyLMM using `pip`, see https://github.com/nickFurlotte/pylmm for details.
The mousegwas will also download METASOFT and run it on the output if there is more than one phenotype.

As part of the data processing, mousegwas can select a subset of the individuals, restricting the number of mice in each strain x sex group or use the average phenotype of all the individuals in such a group. This is controlled by the `-d` option with 0 for average or any other integer for number restriction.

By default LOCO will be used, use the `--noloco` argument to disable it.

A quantile-quantile normalizatin of each phenotype meausrement could be done using the `--qqnorm` argument. 
Other parameters will control the final Manhattan plot, it is a bit unnecessary since the `postprocess_GWAS.R` script will generate more and publication ready figures. 

Nextflow pipeline
-----------------
To execute the scripts in an easy way we included a nextflow pipeline that runs the initial GWAS, the shuffled executions, 
determine a p-value and run the postprocess.
To run coat color phenotype GWAS you can simply install `nextflow`, make sure that `singularity` is installed and run:
```bash
nextflow run TheJacksonLaboratory/mousegwas \
  --yaml https://raw.githubusercontent.com/TheJacksonLaboratory/mousegwas/master/example/coat_color_MDA.yaml \
  --shufyaml https://raw.githubusercontent.com/TheJacksonLaboratory/mousegwas/master/example/coat_color_MDA.yaml \
  --addgwas="--coat_phenotype" --addpostp="--coat_phenotype --colorgroup" --pvalue 0.1 --clusters 1 --outdir coatout -profile singularity,slurm
```
`slurm` can be changed to `pbs` or ignored for local execution.

To regenerate the results in the paper: https://www.biorxiv.org/content/10.1101/2020.10.08.331017v1 :
```bash
nextflow run TheJacksonLaboratory/mousegwas \
  --yaml https://raw.githubusercontent.com/TheJacksonLaboratory/mousegwas/master/example/grooming_nowild.yaml \
  --shufyaml https://raw.githubusercontent.com/TheJacksonLaboratory/mousegwas/master/example/grooming_shuffle.yaml \
  --input https://raw.githubusercontent.com/TheJacksonLaboratory/mousegwas/master/example/grooming_paper_strain_survey_2019_11_21.csv \ \
  --outdir grooming_output --addpostp="--loddrop 0" -profile slurm,singularity
```

