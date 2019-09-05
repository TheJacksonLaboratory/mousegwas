#!/usr/bin/env Rscript
#
# Read csv file, prepare for GEMMA, run GEMMA and plot the results
#
library(gemmawrapper)
library(argparse)
library(yaml)
library(dplyr)
library(readr)
library(R.utils)
library(data.table)
library(snpStats)
library(ggplot2)

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--downsample", '-d', default=0, type="integer",
                    help="Downsample strains to have at most this number of representatives. If 0 (default) average each strain")
parser$add_argument("-i", "--input",
                    help="Input csv file, assume column names. Definitions for which columns to use and how are in the yaml file")
parser$add_argument("-y", "--yaml",
                    help="Yaml file which will include: strain- Name of strain column, phenotypes-phenotypes column names to test, covar-columns to be used as covariates, F1- a dictionary mapping F1 names to strain names, translate- translate strain names to names in the genotypes file")
parser$add_argument("--gemma",
                    help="Gemma executable")
parser$add_argument("-g", "--genotypes", nargs = '+',
                    help="The genotypes csv files as downloaded from phenome.jax.org")
parser$add_argument("--basedir", default=".",
                    help="output directory. Will overwrite existing files")
parser$add_argument("--genes", default=NULL,
                    help="a tab delimited table with SNP ID in first column and gene name in second. Might not include all SNPs")
parser$add_argument("--pylmm", "-p", default="pylmmGWAS.py",
                    help="pylmmGWAS.py executable")
parser$add_argument("--pylmmkinship", "-k", default="pylmmKinship.py",
                    help="pylmmKinship.py executable")
parser$add_argument("--method", "-m", default="GEMMA",
                    help="Which GWAS software to use, possible values: GEMMA|pyLMM")
parser$add_argument("--noloco", default=FALSE, action="store_true",
                    help="Do not use LOCO")
args <- parser$parse_args()

# Load the yaml
yamin <- yaml.load_file(args$yaml)

# Make the output directory
dir.create(args$basedir, recursive = TRUE)
# Read the input table
complete_table <- read_csv(args$input, col_names=TRUE)

# Generate a list of phenotypes
pheno_names <- intersect(yamin$phenotypes, names(complete_table))
# A list of covariates
covar_names <- setdiff(intersect(yamin$covar, names(complete_table)), pheno_names)

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

longfile <- tempfile()
complete.geno <- NULL
for (f in args$genotypes){
  geno <- fread(f)
  geno[, c("major", "minor") := tstrsplit(observed, "/", fixed=TRUE, keep=1:2)]

  if (is.null(complete.geno)){
    if (length(intersect(names(geno), valid_strains))>0){
      addnames <- c(intersect(names(geno), valid_strains), "chr", "bp38", "rs", "major", "minor")
      complete.geno <- geno[,..addnames]
    }
  }else{
    addnames <- c(intersect(names(geno), valid_strains), "chr", "bp38", "rs", "major", "minor")
    geno <- geno[,..addnames]
    setkey(geno, rs)
    setkey(complete.geno, rs)
    complete.geno <- merge(complete.geno, geno, all=TRUE, by=c("chr", "bp38", "rs", "major", "minor"))

  }
  valid_strains <- setdiff(valid_strains, names(complete.geno))
}
complete.geno[, chr:=as.character(chr)]
srdata <- complete.geno[, .(rs, major, minor)]
numeric.geno <- complete.geno[, .("chr", "bp38", "rs", "major", "minor")]
for (cn in setdiff(names(complete.geno), c("chr", "bp38", "rs", "major", "minor"))){
  complete.geno[get(cn)=='H',c(cn) := 1]
  complete.geno[get(cn)==major, c(cn) := 0]
  complete.geno[get(cn)==minor, c(cn) := 2]
  complete.geno[,c(cn) := as.numeric(get(cn))]#as.numeric(ifelse(..cn=='H', 1, ifelse(..cn==major, 0, 2)))]
}

#fwrite(complete.geno, "export_complete_genotypes.csv")
#fwrite(complete.geno[,.(rs, bp38, chr)], "export_snp_data.csv", col.names = FALSE, na="NA")
# Compute the specific strains genomes
strains_genomes <- srdata
sorder = c()
phenos = data.table()
for (p in pheno_names) phenos <- phenos[, eval(p) := numeric()]
covars = data.table()
for (c in covar_names) covars <- covars[, eval(c) := numeric()]
covars <- covars[, isWild := numeric()]

#write_csv(strains, paste0(args$basedir, "/export_strains.csv"))
for (comrow in 1:dim(complete_table)[1]){
  sname <- as.character(complete_table[comrow, yamin$strain])
  rnum <- which(strains$input_name==sname)
  if (length(rnum) == 0){
    print(paste0("Can't find strain data ", sname, " row ",comrow))
    next
  }
  p1n <- strains$p1[rnum]
  p2n <- strains$p2[rnum]
  if (p1n %in% names(complete.geno) & p2n %in% names(complete.geno) & (!p1n %in% yamin$wild)){
    sorder <- c(sorder, sname)
    strains_genomes[, eval(paste0('X',comrow)):=(complete.geno[,..p1n] + complete.geno[,..p2n])/2]
    # Add the phenotypes to the table
    phenos <- rbind(phenos, complete_table[comrow, pheno_names])
    # Add the covariates to the table
    covars <- rbind(covars, cbind(complete_table[comrow, covar_names], tibble(isWild=as.numeric(p1n %in% yamin$wild | p2n %in% yamin$wild))))
  }else{
    print(paste0("Can't find ", p1n," or ", p2n))
  }
}
# Remove SNPs with more than 5% missing data
strains_genomes <- strains_genomes[rowSums(is.na(strains_genomes))<(ncol(strains_genomes)-3)/20,]
# Add 1 to covariates matrix
covars <- model.matrix(as.formula(paste0("~", do.call(paste, c(as.list(covar_names), sep="+")))), covars) #cbind(1, covars)

# Export order of strains used
write.csv(sorder, paste0(args$basedir, "/export_strains_order.csv"))

# Take the betas of each strain and use it to run GEMMA
b <- average_strain(strains_genomes, phenos, covars, args$downsample)

# Run gemma using the helper function with loco

if (args$method == "GEMMA"){
  results_file <- execute_lmm(data.table(b$genotypes), data.table(b$phenotypes),
                              as.data.table(complete.geno[,.(rs, bp38, chr)]),
                              NULL, args$basedir, yamin$eigens, loco=!args$noloco, single=is.null(yamin$eigens) || (yamin$eigens==0))
}else if (args$method == "pyLMM"){
  results_file <- run_pylmm(data.table(b$genotypes), data.table(b$phenotypes),
                                as.data.table(complete.geno[,.(rs, bp38, chr)]),
                                NULL, args$basedir, args$pylmm, args$pylmmkinship, loco=!args$noloco)
}


p <- plot_gemma_lmm(results_file, metasoft=is.null(yamin$eigens) || (yamin$eigens==0), annotations=paste0(args$basedir, "/annotations.csv"))
ggsave(paste0(args$basedir, "/manhattan_plot_p_lrt.pdf"), plot=p, device="pdf", width=16, height=8, units="in")
