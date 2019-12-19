#!/usr/bin/env Rscript
#
# Plot figures and other post-processing
#
# Get the results directory of a run_GEMMA.R script and prepare publication-ready figures (yeah!)
# The script was built and tested on GEMMA only with multiple phenotypes.
# Input also includes the number of clusters to use for clustering. Colors are preset

# Load relevant libraries:
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
#library(ggbiplot)
library(RColorBrewer)
library(viridis)
library(gplots)
library(argparse)
library(gemmawrapper)

parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--outdir", "-o",
                    help="run_GEMMA.R output dir")
parser$add_argument("--plotdir", "-p", default=".",
                    help="Where to put the plots")
parser$add_argument("--clusters", '-c', default=5, type="integer",
                    help="Number of peaks clusters")
args <- parser$parse_args()

# Step 1: Read the color pallete
# Cluster colors
ccols <- brewer.pal(args$clusters, "Dark2")
# Heatmap plot for m-values
hmcol <- viridis(256, option="cividis")

# Read the data
# Read the METASOFT results file. The names of the columns are taken from the phenotypes_order file
phenos <- as.character(read.csv(paste0(args$outdir, "/phenotypes_order.txt"), header = FALSE, skip=1)$V1)
print(phenos)
cnames <- c("rs", "STUDYNUM", "PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE",
            "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q",
            "TAU_SQUARE", paste0(phenos, "_PV"), phenos, "empty")
allgwas <- read_delim(paste0(args$outdir, "/output/all_lmm_associations.assoc.txt"), "\t", col_names = cnames, skip=1, guess_max = 10000)
anno <- read_delim(paste0(args$outdir, "/annotations.csv"), ",", col_names = c("rs", "ps", "chr"), guess_max = Inf)
allgwas <- left_join(allgwas, anno, by="rs") %>% arrange(chr, ps)
# Read the plotting data
#p <- load(paste0(args$outdir, "gwas_object_output.Rdata"))
# Read the genotypes
geno <- as.matrix(read_csv(paste0(args$outdir, "/strains_genotypes_all.csv"), col_types = cols(
  .default = col_double(),
  chr = col_character(),
  rs = col_character(),
  major = col_character(),
  minor = col_character()
)) %>% column_to_rownames(var = "rs") %>% dplyr::select(-chr, -major, -minor))

# We're all set
dir.create(args$plotdir, recursive = TRUE)
# Plot the combined Manhattan plot
p <- plot_gemma_lmm(paste0(args$outdir, "/output/all_lmm_associations.assoc.txt"),
                    annotations = paste0(args$outdir, "/annotations.csv"), metasoft = TRUE,
                    name = "Chromosome",genotypes = geno, namethr = 7, redthr = 7)
ggsave(filename = paste0(args$plotdir, "/Manhattan_plot_all_phenotypes.pdf"),
       plot=p$plot + theme(text=element_text(size=10, family="Times")), device="pdf", dpi="print", width=7.25, height=3.6, units="in")

# Plot each phenotype's Manhattan plot
for (i in 1:length(phenos)){
  pp <- plot_gemma_lmm(Sys.glob(paste0(args$outdir, "/output/lmm_*_pheno_", i, ".assoc.txt")),
                       name = "Chromosome",genotypes = geno, namethr = 7, redthr = 7)
  ggsave(filename = paste0(args$plotdir, "/Manhattan_plot_phenotype_", i, "_", phenos[i], ".pdf"),
         plot=pp$plot + theme(text=element_text(size=10, family="Times")), device="pdf", dpi="print", width=7.25, height=3.6, units="in")
}


