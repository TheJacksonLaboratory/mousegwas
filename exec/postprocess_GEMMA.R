#!/usr/bin/env Rscript
#
# Plot figures and other post-processing
#
# Get the results directory of a run_GEMMA.R script and prepare publication-ready figures (yeah!)
# The script was built and tested on GEMMA only with multiple phenotypes.
# Input also includes the number of clusters to use for clustering. Colors are preset

# Load relevant libraries:
library(plyr)
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
fullw <- 7.25
halfw <- 3.54
height <- 3.54
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

PVE <- read_csv(paste0(args$outdir, "/PVE_GEMMA_estimates.txt"))
# We're all set
dir.create(args$plotdir, recursive = TRUE)
# Plot the combined Manhattan plot
p <- plot_gemma_lmm(paste0(args$outdir, "/output/all_lmm_associations.assoc.txt"),
                    annotations = paste0(args$outdir, "/annotations.csv"), metasoft = TRUE,
                    name = "Chromosome",genotypes = geno, namethr = 15, redthr = 10, maxdist=10000000)
ggsave(filename = paste0(args$plotdir, "/Manhattan_plot_all_phenotypes.pdf"),
       plot=p$plot + theme(text=element_text(size=10, family="Times")), device="pdf", dpi="print",
       width=fullw, height=height, units="in")

# Plot each phenotype's Manhattan plot
for (i in 1:length(phenos)){
  pp <- plot_gemma_lmm(Sys.glob(paste0(args$outdir, "/output/lmm_*_pheno_", i, ".assoc.txt")),
                       name = "Chromosome",genotypes = geno, namethr = 7, redthr = 7, maxdist=10000000)
  ggsave(filename = paste0(args$plotdir, "/Manhattan_plot_phenotype_", i, "_", phenos[i], ".pdf"),
         plot=pp$plot + theme(text=element_text(size=10, family="Times")), device="pdf", dpi="print",
         width=fullw, height=height, units="in")
}

# Cluster the peaks using the m-values
set.seed(49)
pgwas <- allgwas %>% filter(rs %in% p$gwas$rs[p$gwas$ispeak]) %>% column_to_rownames(var = "rs")
pgwas <- as.matrix(pgwas[, phenos])
pcvals <- prcomp(pgwas)
pcmvals <- cbind(pgwas, pcvals$x)
pcperc <- pcvals$sdev^2/sum(pcvals$sdev^2)
kk <- kmeans(pgwas, args$clusters, nstart=5)
pcmvals <- cbind(pcmvals, kk$cluster)
# plot the PCA
bip <- ggbiplot::ggbiplot(pcvals, groups=as.factor(kk$cluster)) + scale_color_manual(name = 'cluster', values=ccols)
ggsave(filename = paste0(args$plotdir, "/PCA_plot.pdf"),
       plot = bip + theme(text=element_text(size=10, family="Times")),
       device="pdf", dpi="print", width=halfw, height=height, units="in")
# Plot the m-value heatmap
clustcol <- tibble(cluster=1:args$clusters, color=ccols)
colrow <- tibble(rs = colnames(pgwas), cluster=kk$clusters) %>% left_join(clustcol, by="cluster") %>% column_to_rownames(var = "rs") %>% dplyr::select(color)
pdf(paste0(args$plotdir, "/all_peaks_heatmap.pdf"), width = fullw, height = height, family = "Times")
heatmap.2(pgwas, col = hmcol,
          Rowv = T, Colv = T, dendrogram = "both", scale="none", trace="none",
          RowSideColors = colrow[,1,drop=T], labRow = NA)
dev.off()

# Plot the PVE estimates with SE
pvep <- ggplot(PVE, aes(reorder(phenotype, -PVE), PVE)) + geom_bar(color="black", fill = RColorBrewer::brewer.pal(3,"Set1")[2],
                                            stat="identity") +
  geom_errorbar(aes(ymin=PVE-PVESE, ymax=PVE+PVESE), width=.2) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(text=element_text(size=10, family="Times"))
ggsave(paste0(args$plotdir, "/PVE_plot.pdf"), plot = pvep, device = "pdf", dpi = "print",
       width = fullw, height = height, units = "in")
