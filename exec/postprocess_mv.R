#!/usr/bin/env Rscript
#
# Plot figures and other post-processing
#
# Get the results directory of a run_GEMMA.R script and prepare publication-ready figures
# The script was built and tested on GEMMA only with multiple phenotypes.
# Input also includes the number of clusters to use for clustering. Colors are preset

# Load relevant libraries:
library(plyr)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gplots)
library(argparse)
library(enrichR)
library(mousegwas)

parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--outdir", "-o",
                    help="run_GEMMA.R output dir")
parser$add_argument("--plotdir", "-p", default=".",
                    help="Where to put the plots")
parser$add_argument("--clusters", '-c', default=7, type="integer",
                    help="Number of peaks clusters")
parser$add_argument("--rotation", "-r", default="",
                    help="comma separated list of phenotypes to plot in the ggbiplot. If empty plots all")
parser$add_argument("--sample", "-s", type="integer", default=10000,
                    help="Number of SNPs to sample for the LD plotting")
parser$add_argument("--names", "-n",
                    help="Translation of the phenotypes to paper names. The csv file should have the columns Group, OriginalName, PaperName")
parser$add_argument("--pvalthr", default=5, type="double",
                    help="p-value threshold for plotting and getting gene lists")
args <- parser$parse_args()

# Step 1: Read the color pallete
# Cluster colors
ccols <- brewer.pal(args$clusters, "Dark2")[1:args$clusters]
# Heatmap plot for m-values
pigr <- RColorBrewer::brewer.pal(name = "PiYG", n=11)
hmcol <- viridis(128)#colorRampPalette(pigr[c(2,5,10)])(128)#viridis(128, option="cividis")
grpcol <- RColorBrewer::brewer.pal(8,"Accent")
fullw <- 7.25
halfw <- 3.54
height <- 3.54
ffam <- "Arial"
# Read the data
# Phenotypes names
phenos <- as.character(read.csv(paste0(args$outdir, "/phenotypes_order.txt"), header = FALSE, skip=1)$V1)
pnames <- read.csv(args$names,row.names = 2)
phenos <- as.character(pnames[phenos, "PaperName", drop=T])
# SNPs annotations
anno <- read_delim(paste0(args$outdir, "/annotations.csv"), ",", col_names = c("rs", "ps", "chr"), guess_max = Inf)
# Read the p-values from all the pasted phenotypes and mv LOCO files
pvalmat <- NULL

# Read the genotypes
geno_t <- read_csv(paste0(args$outdir, "/strains_genotypes_all.csv"), col_types = cols(
  .default = col_double(),
  chr = col_character(),
  rs = col_character(),
  major = col_character(),
  minor = col_character()
))

geno <- as.matrix(geno_t %>% column_to_rownames(var = "rs") %>% dplyr::select(-chr, -bp38, -major, -minor))

PVE <- read_csv(paste0(args$outdir, "/PVE_GEMMA_estimates.txt"))
PVE <- left_join(PVE, as_tibble(pnames, rownames="phenotype"), by = ("phenotype"))

# We're all set
dir.create(args$plotdir, recursive = TRUE)
set.seed(490)

# Plot each phenotype's Manhattan plot
lilp <- vector("list", length(phenos))
allpeaks <- NULL
allsnps <- NULL
all_ispeak <- NULL
all_choose <- NULL
for (i in 1:length(phenos)){
  pp <- plot_gemma_lmm(paste0(args$outdir, "/output/lmm_pheno_", i, "_all_LOCO.assoc.txt"),
                       name = "Chromosome",genotypes = geno, namethr = args$pvalthr, redthr = args$pvalthr, maxdist=10000000,
                       corrthr=0.4)
  pname <- phenos[i]
  lilp[[pname]] <- pp
  allsnps <- pp$gwas$rs
  if (is.null(pvalmat)){
    pvalmat <- pp$gwas %>% mutate(!!(pname):=P) %>% dplyr::select(rs, !!(pname))
    all_ispeak <- pp$gwas %>% mutate(!!(pname):=ispeak) %>% dplyr::select(rs, !!(pname))
    all_choose <- pp$gwas %>% mutate(!!(pname):=choose) %>% dplyr::select(rs, !!(pname))
  }else{
    pvalmat <- left_join(pvalmat, pp$gwas %>% mutate(!!(pname):=P) %>% dplyr::select(rs, !!(pname)), by="rs")
    all_ispeak <- left_join(all_ispeak, pp$gwas %>% mutate(!!(pname):=ispeak) %>% dplyr::select(rs, !!(pname)), by="rs")
    all_choose <- left_join(all_choose, pp$gwas %>% mutate(!!(pname):=choose) %>% dplyr::select(rs, !!(pname)), by="rs")
    print(head(all_choose))
  }
  allpeaks <- c(allpeaks, pp$gwas$rs[pp$gwas$ispeak])
  ggsave(filename = paste0(args$plotdir, "/Manhattan_plot_phenotype_", i, "_", phenos[i], ".pdf"),
         plot=pp$plot + theme(text=element_text(size=10, family=ffam)), dpi="print", device = cairo_pdf,
         width=fullw, height=height, units="in")
}
# Plot the groups MAnhattan plots
grouplist <- c()
for (grpf in Sys.glob(paste0(args$outdir,"/output/lmm_phenotypes_*_all_LOCO.assoc.txt"))){
  pp <- plot_gemma_lmm(grpf, name = "Chromosome",genotypes = geno, namethr = args$pvalthr, redthr = args$pvalthr,
                       maxdist=10000000, corrthr=0.4)
  pname <- gsub(".*phenotypes_(.*)_all_LOCO.assoc.txt", "\\1", grpf)
  grouplist <- c(grouplist, pname)
  allres <- read_delim(grpf, "\t", guess_max = 1000000) %>% mutate(!!(pname) := p_score) %>% dplyr::select(rs, !!(pname))
  pvalmat <- left_join(pvalmat, pp$gwas %>% mutate(!!(pname) := P) %>% dplyr::select(rs, !!(pname)), by="rs")
  all_ispeak <- left_join(all_ispeak, pp$gwas %>% mutate(!!(pname):=ispeak) %>% dplyr::select(rs, !!(pname)), by="rs")
  all_choose <- left_join(all_choose, pp$gwas %>% mutate(!!(pname):=choose) %>% dplyr::select(rs, !!(pname)), by="rs")

  allpeaks <- c(allpeaks, pp$gwas$rs[pp$gwas$ispeak])
  ggsave(filename = paste0(args$plotdir, "/Manhattan_plot_phenotypes_", pname, ".pdf"),
         plot=pp$plot + theme(text=element_text(size=10, family=ffam)), dpi="print", device = cairo_pdf,
         width=fullw, height=height, units="in")
}

# Cluster the peaks using the P values

pgwas <- pvalmat %>% filter(rs %in% allpeaks) %>% column_to_rownames(var = "rs")
print(sum(is.na(pgwas)))
pgwas[is.na(pgwas)] = 0
pgwas <- as.matrix(pgwas)
#pcvals <- prcomp(pgwas)
#if (args$rotation != ""){
#  pcvals$rotation <- pcvals$rotation[strsplit(args$rotation, ",")[[1]],,drop=F]
#}
#pcperc <- pcvals$sdev^2/sum(pcvals$sdev^2)
kk <- kmeans(pgwas, args$clusters)
# plot the PCA
#bip <- ggbiplot::ggbiplot(pcvals, groups=as.factor(kk$cluster)) + scale_color_manual(name = 'cluster', values=ccols) + theme_bw() + theme(legend.position = "none")

#ggsave(filename = paste0(args$plotdir, "/PCA_plot.pdf"),
#       plot = bip + theme(text=element_text(size=10, family=ffam)),
#       device = cairo_pdf, dpi="print", width = halfw, height = height, units="in")


# Plot the m-value heatmap
# Order the columns
rowarr <- NULL
for (i in 1:args$cluster){
  cc <- hclust(dist(pgwas[kk$cluster==i,, drop=F]))
  rowarr <- c(rowarr, row.names(pgwas)[kk$cluster==i][cc$order])
}
pgwas <- pgwas[rowarr,]
clustcol <- tibble(cluster=1:args$clusters, color=ccols)
colrow <- tibble(rs = rownames(pgwas), cluster=kk$cluster[rowarr]) %>% left_join(clustcol, by="cluster") %>% column_to_rownames(var = "rs") %>% dplyr::select(color)
grptocol <- tibble(Group=unique(PVE$Group), color=grpcol[1:length(unique(PVE$Group))])
colcol <- PVE %>% left_join(grptocol) %>% column_to_rownames(var="PaperName") %>% dplyr::select(color)
colcol <- colcol$color
# Add the groups to the colors
for (g in grouplist){
  colcol <- c(colcol, grptocol[grepl(g, grptocol$Group),]$color)
}
cairo_pdf(paste0(args$plotdir, "/all_peaks_heatmap.pdf"), width = fullw, height = height+1, family = ffam)
heatmap.2(pgwas, col = hmcol,
          Rowv = F, Colv = T, dendrogram = "col", scale="none", trace="none",
          RowSideColors = colrow[,1,drop=T], ColSideColors = colcol, labRow = NA,
          hclustfun = function(x) hclust(x, method="average"),
          margins=c(12,8),srtCol=45, key=T, density.info = "none")
dev.off()


# Plot the PVE estimates with SE
pvep <- ggplot(PVE, aes(reorder(PaperName, -PVE), PVE, fill=Group)) + geom_bar(color="black", stat="identity") +
  scale_fill_manual(values = grpcol) +
  geom_errorbar(aes(ymin=PVE-PVESE, ymax=PVE+PVESE), width=.2) +
  xlab("Phenotype") +
  theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(text=element_text(size=10, family=ffam))
ggsave(paste0(args$plotdir, "/PVE_plot.pdf"), plot = pvep, device = cairo_pdf, dpi = "print",
       width = fullw, height = height, units = "in")


# Plot the metasoft manhattan plot with clusters colors
# Add the cluster number to the pwas object
#p$pwas <- p$pwas %>% left_join(tibble(rs = rownames(pgwas), cluster=as.factor(kk$cluster)), by="rs")
#p$pwas <- p$pwas %>% dplyr::select(-cluster) %>% left_join(filter(p$pwas,ispeak)%>%dplyr::select(choose, cluster),by="choose")

# Recolor the second layer with the clusters colors
#pnoname <- p$plot
#pnoname$layers <- pnoname$layers[1:2]
#nolab <- p$plot
#nolab$layers <- nolab$layers[1]
#ggsave(filename = paste0(args$plotdir, "/replot_Manhattan_clusters_all.pdf"),
#       plot = pnoname + ggnewscale::new_scale_color() +
#         geom_point(aes(color=p$pwas$cluster), size=0.9) +
#         scale_color_manual(values=ccols) +
#         ggnewscale::new_scale_color() +
#         geom_point(aes(alpha=p$pwas$ispeak), size=1.2, color="black") +
#         scale_alpha_manual(values = c(0,1)) +
#         ggnewscale::new_scale_color() +
#         geom_point(aes(color=p$pwas$cluster, alpha=p$pwas$ispeak), size=0.9) +
#         scale_color_manual(values=ccols) +
#         scale_alpha_manual(values = c(0,1)) +
#         theme(text=element_text(size=10, family=ffam)),
#       device=cairo_pdf, dpi="print", width=fullw, height=height, units="in")
# Plot each cluster's Manhattanplot
#p$pwas <- p$pwas %>% select(-cluster) %>% left_join(filter(p$pwas,ispeak)%>%select(choose, cluster),by="choose")
#for (k in 1:args$clusters){
#  ggsave(filename = paste0(args$plotdir, "/replot_Manhattan_cluster_", k, ".pdf"),
#         plot = nolab %+% p$pwas[p$pwas$cluster==k | is.na(p$pwas$cluster),] + ggnewscale::new_scale_color() +
#           geom_point(aes(color=cluster), size=0.9) +
#           scale_color_manual(values=ccols[k]) +
#           ggnewscale::new_scale_color() +
#           geom_point(aes(alpha=ispeak), size=1.2, color="black") +
#           scale_alpha_manual(values = c(0,1)) +
#           ggnewscale::new_scale_color() +
#           geom_point(aes(color=cluster, alpha=ispeak), size=0.9) +
#           scale_color_manual(values=ccols[k]) +
#           scale_alpha_manual(values = c(0,1)) +
#           theme(text=element_text(size=10, family=ffam)),
#         device=cairo_pdf, dpi="print", width=fullw, height=height, units="in")
#}


# Plot the LD drop figure
comp_LD_2 <- function(genotypes, c, maxdist = 2500000, MAF=0.1, miss=0.1){
  gen <- genotypes %>% filter(chr==c, !is.na(rs)) %>% arrange(bp38)
  gen <- gen[(rowSums(gen[,6:ncol(gen)]==0, na.rm=T) > (dim(gen)[2]-5)*MAF &
                rowSums(gen[,6:ncol(gen)]>0, na.rm=T) > (dim(gen)[2]-5)*MAF &
                rowSums(is.na(gen[,6:ncol(gen)])) < (dim(gen)[2]-5)*miss),]
  allcor = NULL
  for (i in 1:(nrow(gen)-1)){
    tmat <- as.data.frame(base::t(gen %>% filter(between(bp38, gen$bp38[i], gen$bp38[i]+maxdist)) %>%
                                    dplyr::select(-bp38, -chr, -major, -minor) %>% column_to_rownames(var = "rs")))
    if (ncol(tmat)<2) next()
    ct <- cor(tmat[,1], tmat[, 2:ncol(tmat), drop=F], use = 'pairwise.complete.obs', method="pearson")
    allcor <- rbind(allcor, as_tibble(t(ct), rownames = "SNP2") %>% mutate(r_sq = V1^2, SNP1=gen$rs[i]) %>% dplyr::select(SNP1, SNP2, r_sq))
  }
  return(allcor)
}
geno_s <- geno_t %>% filter(chr != "X", chr != "Y", chr != "MT") %>% sample_n(args$sample)
allchr = NULL
for (chr in unique(geno_s$chr)){
  allchr <- rbind(allchr, comp_LD_2(geno_s, chr))
}
allchr <- allchr %>% left_join(dplyr::select(geno_s, rs, bp38), by=c("SNP1" = "rs")) %>%
  left_join(dplyr::select(geno_s, rs, bp38), by=c("SNP2" = "rs")) %>% mutate(dist = bp38.y-bp38.x)
avgwin = 5000
allchr$distc <- (cut(allchr$dist, breaks=seq(from=min(allchr$dist)-1,to=max(allchr$dist)+1,by=avgwin)))
allavg <- allchr %>% group_by(distc) %>% summarise(avdist=mean(dist),avr_sq=mean(r_sq, na.rm = T)) %>% ungroup()
pld <- ggplot(allavg, aes(avdist/1000000, avr_sq)) + geom_smooth(color=RColorBrewer::brewer.pal(3,"Set1")[3], method="loess", span=0.3, se=FALSE)+
  xlim(c(0, 2.5)) + labs(x="Distance (Mbp)",y=expression("Average LD"~(r^{2})))
ggsave(filename = paste0(args$plotdir, "/plot_LD_drop.pdf"), plot=pld + theme_bw() + theme(
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  text=element_text(size=10, family=ffam)
),
device=cairo_pdf, dpi="print", width=halfw, height=height, units="in"
)


# Plot MAF histogram
#mafdat <- tibble(rs = geno_t$rs, maf = rowSums(geno_t[,-1:-5])/(2*(ncol(geno_t)-5)))
#mafdat$maf <- pmin(mafdat$maf, 1-mafdat$maf)
#mafdat <- mafdat %>% filter(rs %in% allsnps)
#mafp <- ggplot(mafdat, aes(maf, fill=choose==0, color=choose==0)) + geom_histogram(binwidth = 1/(ncol(geno_t)-5)) + xlim(c(0,0.5)) +
#  scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")[3:4], name="", labels=c("All","Peak")) +
#  scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Paired")[3:4], name="", labels=c("All","Peak")) +
#  xlab("MAF") +
#  theme_bw() + theme(legend.position=c(0.15,0.9),
#                     panel.border = element_blank(),
#                     panel.grid.major.x = element_blank(),
#                     panel.grid.minor.x = element_blank(),
#                     panel.grid.major.y = element_blank(),
#                     panel.grid.minor.y = element_blank(),
#                     text=element_text(size=10, family=ffam)
#  )
#ggsave(filename = paste0(args$plotdir, "/plot_MAF_hist.pdf"), plot=mafp,
#       device=cairo_pdf, dpi="print", width=halfw, height=height, units="in")



# Get the genes related to each cluster, print them and run enrichR
ext_peak <- function(snps, all_i, all_c, maxdist=2000000){
  csum <- snps %>% filter(rs %in% all_i$rs[rowSums(all_i %>% select(-rs))>0])
  csum$minps <- csum$ps
  csum$maxps <- csum$ps
  for (c in names(all_c %>% select(-rs))){
    tmps <- left_join(snps, select(all_c, rs, !!(c)), by="rs") %>% filter(!!(c) > 0)
    tmpc <- tmps %>% group_by_at(c) %>% dplyr::summarize(maxps = max(ps), minps = min(ps)) %>% ungroup() %>% left_join(tmps, by=c) %>% select(rs, minps, maxps)
    csum <- csum %>% left_join(tmpc, by="rs", suffix=c("",".x")) %>% mutate(maxps = pmin(maxps, maxps.x), minps = pmax(minps, minps.x))
  }
  csum %>% mutate(maxps = pmin(maxps, ps + maxdist), minps = pmax(maxps, ps - maxdist))

}
#pg <- ext_peak(left_join(p$gwas, dplyr::select(p$pwas, rs, cluster),by="rs"))
#pg <- pg %>% left_join(dplyr::select(allgwas, rs, STAT1_RE2, STAT2_RE2))
#ggsave(filename = paste0(args$plotdir, "/peak_width_hist.pdf"),
#       plot = ggplot(filter(pg, ispeak), aes(log10(maxps-minps))) + geom_histogram(bins = 20) +
#         theme_bw() + theme(text=element_text(size=10, family=ffam)),
#       device=cairo_pdf, dpi="print", width=halfw, height=height, units="in")
pg <- anno %>% filter(rs %in% allsnps) %>% join(tibble(rs = rownames(pgwas), cluster=kk$cluster[rowarr]))
pg <- ext_peak(pg, all_ispeak, all_choose)
allgenes = NULL
dbs <- listEnrichrDbs()

for (k in 1:args$clusters){
  affgen <- get_genes(pg[pg$ispeak==T & pg$cluster==k,], dist=1000)
  allgenes <- rbind(allgenes, affgen)
  write_csv(affgen, path = paste0(args$plotdir, "/genes_for_cluster_", k, ".csv"))
  # Run enrichr
  enrr <- enrichr(unique(affgen$mgi_symbol[!(grepl(pattern = "^Gm", x =  affgen$mgi_symbol) | grepl("Rik$", affgen$mgi_symbol) | affgen$mgi_symbol=="")]), dbs$libraryName)
  for (d in dbs$libraryName){
    if (dim(enrr[[d]])[2] > 1){
      rtb <- as_tibble(enrr[[d]]) %>% filter(Adjusted.P.value <= 0.05 )
      if (nrow(rtb) > 0)
        write_csv((rtb %>% mutate(total.genes= length(affgen$mgi_symbol[!(grepl(pattern = "^Gm", x =  affgen$mgi_symbol) | grepl("Rik$", affgen$mgi_symbol) | affgen$mgi_symbol=="")]), cluster=k, library=d)),
                  path = paste0(args$plotdir, "/enrichR_cluster_", k, "_db_", d, "p005.csv"))

    }
  }
}

# Do this for each phenotype
for (i in 1:length(lilp)){
  pp <- lilp[[i]]
  expp <- ext_peak(pp$gwas)
  affgen <- get_genes(expp[expp$ispeak==T,], dist=1000)
  if (nrow(affgen) > 0){
    write_csv(affgen, path = paste0(args$plotdir, "/genes_for_hoenotype_", i, ".csv"))
    # Run enrichr
    enrr <- enrichr(unique(affgen$mgi_symbol[!(grepl(pattern = "^Gm", x =  affgen$mgi_symbol) | grepl("Rik$", affgen$mgi_symbol) | affgen$mgi_symbol=="")]), dbs$libraryName)
    for (d in dbs$libraryName){
      if (length(dim(enrr[[d]]))>1 && dim(enrr[[d]])[2] > 1){
        rtb <- as_tibble(enrr[[d]]) %>% filter(Adjusted.P.value <= 0.05)
        if (nrow(rtb) > 0)
          write_csv((rtb %>% mutate(total.genes= length(affgen$mgi_symbol[!(grepl(pattern = "^Gm", x =  affgen$mgi_symbol) | grepl("Rik$", affgen$mgi_symbol) | affgen$mgi_symbol=="")]), phenotype = pnames$PaperName[i], library=d)),
                    path = paste0(args$plotdir, "/enrichR_phenotype_", i, "_db_", d, "p005.csv"))

      }
    }
  }
}

# Plot markers density
chrord <- c("X", 19:1)
densp <- geno_t %>% filter(chr!="Y", chr!="MT")  %>%
  ggplot(aes(bp38/1000000, factor(chr, levels=chrord))) +
  geom_bin2d(binwidth=1, drop=T) + xlab("Position (Mbp)") + ylab ("Chromosome") +
  scale_fill_viridis(name=expression(frac('markers', '1 Mbp'))) +
  ggnewscale::new_scale_fill() +
  geom_point(data = filter(pg, ispeak, chr!="Y", chr!="MT"),
             aes(ps/1000000, factor(chr, levels=chrord), fill=factor(cluster)),
             show.legend = F, inherit.aes = F, color = "black",
             position = position_jitter(height=0.15, width=0), shape=21) +
  scale_fill_manual(values = ccols) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(size=10, family=ffam))
ggsave(filename = paste0(args$plotdir, "/Chromosome_density_plot.pdf"), plot = densp,
       device=cairo_pdf, dpi="print", width=fullw, height=height, units="in")
