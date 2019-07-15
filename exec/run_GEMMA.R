#!/usr/bin/env Rscript
#
# Read csv file, prepare for GEMMA, run GEMMA and plot the results
#

library(argparse)
library(yaml)
library(dplyr)
library(readr)
library(R.utils)
library(data.table)
library(snpStats)

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input",
                    help="Input csv file, assume column names. Definitions for which columns to use and how are in the yaml file")
parser$add_argument("-y", "--yaml",
                    help="Yaml file which will include: strain- Name of strain column, phenotypes-phenotypes column names to test, covar-columns to be used as covariates, F1- a dictionary mapping F1 names to strain names, translate- translate strain names to names in the genotypes file")
parser$add_argument("--gemma",
                    help="Gemma executable")
parser$add_argument("-g", "--genotypes", nargs = '+',
                    help="The genotypes csv files as downloaded from phenome.jax.org")

args <- parser$parse_args()

# Load the yaml
yamin <- yaml.load_file(args$yaml)
# Read the input table
complete_table <- read_csv(args$input, col_names=TRUE)
# Get the strain names
strains <- complete_table %>% select(input_name=one_of(yamin$strain)) %>% unique() %>% mutate(p1=input_name, p2=input_name)
# Translate F1
a <- tibble(input_name=character(), p1=character(), p2=character())
for (fl in yamin$F1){
  a <- add_row(a, input_name=names(fl), p1=unlist(fl)[1], p2=unlist(fl)[2])
}
# Combine the parents in the strains table
strains <- strains %>% left_join(a, by="input_name", suffix=c(".orig","")) %>%
  mutate(p1=ifelse(is.na(p1), p1.orig, p1), p2=ifelse(is.na(p2), p2.orig, p2)) %>% select(input_name, p1, p2)

# Translate names to genotype names
t <- tibble(name=character(), gen_name=character())
for (tr in yamin$translate){
  t <- add_row(t, name=names(tr), gen_name=unlist(tr)[1])
}
# Replace relevant names in the strains table
strains <- strains %>% left_join(t, by=c("p1" = "name")) %>% mutate(p1=ifelse(is.na(gen_name), p1, gen_name)) %>% select(input_name, p1, p2)
strains <- strains %>% left_join(t, by=c("p2" = "name")) %>% mutate(p2=ifelse(is.na(gen_name), p2, gen_name)) %>% select(input_name, p1, p2)
valid_strains <- unique(c(strains$p1, strains$p2))

# For each genotype file read it and add to the long file, use only genotypes in the input
# Read the genotype csv file
#long_form <- tibble(rs=character(), sample=character(), genotype=integer())
srdata <- NULL
longfile <- tempfile()
complete.geno <- NULL
for (f in args$genotypes){
  geno <- fread(f)
  geno[, c("major", "minor") := tstrsplit(observed, "/", fixed=TRUE, keep=1:2)]
  if (is.null(srdata)){
    srdata <- geno[,.(rs,major, minor)]
  }
  if (is.null(complete.geno)){
    complete.geno <- geno[,.(chr, bp38, rs, major, minor, intersect(names(geno), valid_strains))]
  }else{
    addnames <- intersect(names(geno), valid_strains)
    complete.geno <- cbind(complete.geno, geno[,..addnames])
  }
  valid_strains <- setdiff(valid_strains, names(complete.geno))
  print(dim(complete.geno))
  print(names(complete.geno))
  #long_form <- melt(geno, id.vars=c("rs", "major"), measure.vars=intersect(valid_strains, setdiff(names(geno), c("chr", "bp38", "rs", "major", "minor", "observed", "dbsnp142annot", "requested"))))
  #long_form[,value:=ifelse(value=='H', 2, ifelse(value==major, 1, 3))]
  #long_form[,conf=1.00]
  #fwrite(long_form[,.(rs, variable, value, conf)], longfile, append=TRUE, col.names=FALSE)
}
print(colSums(is.na(complete.geno)))
for (cn in setdiff(names(complete.geno), c("chr", "bp38", "rs", "major", "minor"))){
  complete.geno[,c(cn) := as.numeric(ifelse(..cn=='H', 2, ifelse(..cn==major, 1, 3)))]
}
print(head(complete.geno))
print(summary(complete.geno))


# Compute the specific strains genomes
strains_genomes <- srdata
for (rnum in 1:nrow(strains)){
  p1n <- strains$p1[rnum]
  p2n <- strains$p2[rnum]
  if (p1n %in% names(complete.geno) & p2n %in% names(complete.geno)){
    strains_genomes[, eval(strains$input_name[rnum]):=(complete.geno[,..p1n] + complete.geno[,..p2n])/2]
  }else{
    print(paste0("Can't find ", p1n," or ", p2n))
  }
}
fwrite(strains_genomes, "export_strains_genotypes.csv")
# Arrange the SNPs data
#setnames(srdata, c("chr", "bp38", "rs", "major", "minor"), c("chromosome", "position", "rname", "A1", "A2"))
#srdata <- as.data.frame(srdata)
#rownames(srdata) <- srdata$rname
#srdata <- srdata[,c("chromosome", "position", "A1", "A2")]

# Read the long file using snpStats
#snps <- read.long(longfile, fields=c(snp=1, sample=2, genotype=3, confidence=4),
#                  gcodes=c("1", "2", "3"))

#Impute missing genotypes




