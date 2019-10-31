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
library(ggplot2)
library(GGally)

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
parser$add_argument("--missing", default=0.05, type="double",
                    help="Maximal fraction of missing data for marker")
parser$add_argument("--MAF", default=0.1, type="double",
                    help="Minimal value for minor allele frequency")
parser$add_argument("--header", default="GWAS results",
                    help="Manhattan plot header")
parser$add_argument("--shuffle", default=FALSE, action="store_true",
                    help="Shuffle the phenotypes between the different individuals")
parser$add_argument("--seed", type="integer", default=100,
                    help="If shuffle is true, set the seed to avoid repeating the same results but have the option to rerun")
parser$add_argument("--qqnorm", default=FALSE, action="store_true",
                    help="QQNORM each phenotype before analyzing")
parser$add_argument("--genedist", default=1000000, type="integer",
                    help="gene distance from SNP, for gene reporting")
parser$add_argument("--snpthr", default=7, type="double",
                    help="P threshold for gene reporting")
parser$add_argument("--namethr", default=10, type="double",
                    help="Print gene names above this threshold")
parser$add_argument("--coat_covar", default=FALSE, action="store_true",
                    help="Use coat color as defined in yaml file as a covariate")
parser$add_argument("--coat_phenotype", default=FALSE, action="store_true",
                    help="Use coat color as defined in yaml file as a phenotype")

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

# Read the coat from the yaml file
coat_table <- tibble(strain=character(), coat=character())
if (args$coat_phenotype |args$coat_covar){
  for (ct in yamin$coat){
    coat_table <- add_row(coat_table, strain=names(ct), coat=unlist(ct)[1])
  }
  coat_table_mm <- model.matrix(~coat, coat_table)[,-1]
  row.names(coat_table_mm) <- coat_table$strain
}

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

# Compute the specific strains genomes
strains_genomes <- srdata
sorder = c()
phenos = data.table()
for (p in pheno_names) phenos <- phenos[, eval(p) := numeric()]
covars = data.table()
for (c in covar_names) covars <- covars[, eval(c) := numeric()]
covars <- covars[, isWild := numeric()]
sexvec <- c()

# Keep old not found strains to not repeat error messages
notfounds = c()

# Read the input phenotypic data, write the genotypes and phenotypes in the proper tables
for (comrow in 1:dim(complete_table)[1]){
  sname <- as.character(complete_table[comrow, yamin$strain])
  rnum <- which(strains$input_name==sname)
  if (length(rnum) == 0){
    print(paste0("Can't find strain data ", sname, " row ",comrow))
    next
  }
  p1n <- strains$p1[rnum]
  p2n <- strains$p2[rnum]
  if (p1n %in% names(complete.geno) & p2n %in% names(complete.geno)){
    sorder <- c(sorder, sname)
    strains_genomes[, eval(paste0('X',comrow)):=(complete.geno[,..p1n] + complete.geno[,..p2n])/2]
    # Add the phenotypes to the table

    ct <- if (p1n==p2n) as.character(coat_table %>% filter(strain==p1n) %>% select(coat)) else as.character(coat_table %>% filter(strain==sname) %>% select(coat))
    if (is.null(ct) & (args$coat_phenotype | args$coat_covar)){
      print(paste0("Can't find coat color for ", p1n, " or ", sname))
    }
    prow <- complete_table[comrow, pheno_names]
    if (args$coat_phenotype){
      cname <- if (p1n==p2n) p1n else sname
      prow <- cbind(complete_table[comrow, pheno_names], coat_table_mm[cname,,drop=F])
    }
    phenos <- rbind(phenos, prow, fill=TRUE)
    # Add the covariates to the table
    crow <- cbind(complete_table[comrow, covar_names], tibble(isWild=as.numeric(p1n %in% yamin$wild | p2n %in% yamin$wild)))
    if (args$coat_covar){
      crow <- cbind(crow, coat=ct)
    }
    covars <- rbind(covars, crow, fill=TRUE)
    sexvec <- c(sexvec, complete_table[comrow, yamin$sex])
  }else{
    if (p1n==p2n){
      if (!p1n %in% notfounds){
        print(paste0("Can't find strain ",p1n))
        notfounds <- c(notfounds, p1n)
      }
    }else{
      if (!(p1n %in% notfounds & p2n %in% notfounds)){
        print(paste0("Can't find ", p1n," or ", p2n))
        notfounds <- c(notfounds, p1n, p2n)
      }
    }
  }
}

# Compute the covariate matrix
if (length(covar_names) > 0){
  if (args$coat_covar){
    covar_names <- c(covar_names, "coat")
  }
  covars <- model.matrix(as.formula(paste0("~", do.call(paste, c(as.list(covar_names), sep="+")))), covars)
}else{
  covars = NULL
}

# scale to have mean=0 var=1
phenos <- scale(phenos)

# Remove columns with NaNs
for (c in ncol(phenos):1){ if (all(is.na(phenos[,c]))) phenos <- phenos[,-c, drop=F]}

if (args$shuffle){
  set.seed(args$seed)
  phenos <- phenos[sample(nrow(phenos)),]
}

# Take the betas of each strain and use it to run GEMMA
b <- average_strain(strains_genomes, phenos, covars, args$downsample, sexvec, sorder)

# Print the phenotypes order
write.csv(colnames(b$phenotypes), file=paste0(args$basedir, "/phenotypes_order.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

# Plot correlations between phenotypes

#ppr <- ggpairs(cbind(tibble(Strain=used_strains), tibble(b$phenotypes)))
print(head(b$phenotypes))
ppr <- ggpairs(as.data.frame(b$phenotypes))
ggsave(paste0(args$basedir, "/phenotype_correlations.pdf"), plot=ppr, device="pdf", width=16, height=16, units="in")

# Remove SNPs with more than 5% missing data and 5% MAF
b$genotypes <- b$genotypes[rowSums(is.na(b$genotypes))<=(ncol(b$genotypes)-3)*args$missing &
                             rowSums(b$genotypes==0)>=(ncol(b$genotypes)-3)*args$MAF &
                             rowSums(b$genotypes==2)>=(ncol(b$genotypes)-3)*args$MAF,]

# Normalize the phenotypes
if (args$qqnorm){
  for (r in 1:ncol(b$phenotypes)){
    b$phenotypes[,r] <- qqnorm(as.data.frame(b$phenotypes)[,r], plot=F)$x
  }
}
# Export order of strains used
write.csv(sorder[b$indices], paste0(args$basedir, "/export_strains_order.csv"), quote = FALSE, col.names = FALSE)

# Run gemma/pylmm using the helper function
if (args$method == "GEMMA"){
  results_file <- execute_lmm(data.table(b$genotypes), data.table(b$phenotypes),
                              as.data.table(complete.geno[,.(rs, bp38, chr)]),
                              b$covars, args$basedir, yamin$eigens, loco=!args$noloco, single=is.null(yamin$eigens) || (yamin$eigens==0))
  # Run no LOCO to get the unified heritability for each phenotype
  if (!args$noloco){
    all_res <- execute_lmm(data.table(b$genotypes), data.table(b$phenotypes),
                                as.data.table(complete.geno[,.(rs, bp38, chr)]),
                                b$covars, args$basedir, yamin$eigens, loco=FALSE, single=TRUE)
    # Extract the VPE values for each phenotype
  }
  allVPE = data.table(PVE=numeric(), PVESE=numeric(), Vg=numeric(), Ve=numeric())
  for (n in 1:dim(b$phenotypes)[2]){
    fname <- paste0(args$basedir, "/output/lmm_all_phenotype_", n, ".log.txt")
    if (file.exists(fname)){
      sigs <- get_sigmas(fname)
      allVPE <- rbind(allVPE, sigs)
    }
  }
  fwrite(allVPE, file=paste0(args$basedir, "/PVE_GEMMA_estimates.txt"))
}else if (args$method == "pyLMM"){
  results_file <- run_pylmm(data.table(b$genotypes), data.table(b$phenotypes),
                                as.data.table(complete.geno[,.(rs, bp38, chr)]),
                                b$covars, args$basedir, args$pylmm, args$pylmmkinship, loco=!args$noloco)
}

is.metasoft <- TRUE
if (args$method=="GEMMA"){
  if (!is.null(args$eigens) && args$eigens>0) is.metasoft=FALSE
  if (ncol(b$phenotypes)==1) is.metasoft=FALSE
}else if (args$method=="pyLMM"){
  if (ncol(b$phenotypes)==1) is.metasoft=FALSE
}

# Read the genes file
genes = NULL
if (!is.null(args$genes)){
  genes <- fread(file=args$genes, col.names = c("rs", "gene_name"), header = FALSE)
}

pval_thr <- args$snpthr

# Manhattan plot
p <- plot_gemma_lmm(results_file, genes=genes, name=args$header, metasoft=is.metasoft, pyLMM=args$method=="pyLMM" && ncol(b$phenotypes)==1,
                    annotations=paste0(args$basedir, "/annotations.csv"), namethr=args$namethr, redthr=pval_thr)
ggsave(paste0(args$basedir, "/manhattan_plot_p_lrt.pdf"), plot=p$plot, device="pdf", width=16, height=8, units="in")

# Read the significant SNPs and grab their related genes
#affgen <- get_genes(p$gwas[p$gwas$P>args$snpthr,], dist=args$genedist)
#fwrite(merge(data.table(affgen), data.table(p$gwas), by="rs"), paste0(args$basedir, "/genes_dist_", args$genedist, "_pval_", args$snpthr, ".csv"))
