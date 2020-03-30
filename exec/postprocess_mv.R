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
library(yaml)
library(mousegwas)

parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--outdir", "-o",
                    help = "run_GEMMA.R output dir")
parser$add_argument("--plotdir", "-p", default = ".",
                    help = "Where to put the plots")
parser$add_argument("-y", "--yaml",
                    help="Yaml file which includes the groups under phenotypes")
parser$add_argument("--clusters",
                    '-c',
                    default = 7,
                    type = "integer",
                    help = "Number of peaks clusters")
parser$add_argument("--rotation", "-r", default = "",
                    help = "comma separated list of phenotypes to plot in the ggbiplot. If empty plots all")
parser$add_argument("--sample",
                    "-s",
                    type = "integer",
                    default = 10000,
                    help = "Number of SNPs to sample for the LD plotting")
parser$add_argument("--pvalthr",
                    default = 5,
                    type = "double",
                    help = "p-value threshold for plotting and getting gene lists")
parser$add_argument("--nomv",
                    default = FALSE,
                    action = "store_true",
                    help = "Ignore multivariate results and use only single phenotypes")
parser$add_argument("--inrich", default = "inrich",
                    help = "INRICH exe utable and custom parameters, deafult: inrich")
parser$add_argument("--inrich_i",
                    "-i",
                    type = "double",
                    default = 5,
                    help = "Minimal group size for inrich (-i), default 5")
parser$add_argument("--inrich_j",
                    "-j",
                    type = "double",
                    default = 200,
                    help = "Maximal group size for inrich (-j), default 200")
parser$add_argument("--drop",
                    type = "double",
                    default = 2,
                    help = "Log p-value drop to include in peak")
parser$add_argument("--peakdist",
                    type = "integer",
                    default = 1000000,
                    help = "Half the interval width")
parser$add_argument("--ldpeakdist",
                    type = "integer",
                    default = 10000000,
                    help = "Peak maximal width")
parser$add_argument("--peakcorr",
                    type = "double",
                    default = 0.25,
                    help = "SNPs R^2 correlation cutoff for peak determination")

args <- parser$parse_args()

# Step 1: Read the color pallete
# Cluster colors
ccols <- brewer.pal(args$clusters, "Dark2")[1:args$clusters]
# Heatmap plot for m-values
pigr <- RColorBrewer::brewer.pal(name = "PiYG", n = 11)
hmcol <-
  viridis(128)#colorRampPalette(pigr[c(2,5,10)])(128)#viridis(128, option="cividis")
grpcol <- RColorBrewer::brewer.pal(8, "Accent")
fullw <- 7.25
halfw <- 3.54
height <- 3.54
ffam <- "Arial"
# Read the data
# Phenotypes names
phenos <-
  as.character(read.csv(
    paste0(args$outdir, "/phenotypes_order.txt"),
    header = FALSE,
    skip = 1
  )$V1)
#pnames <- read.csv(args$names, row.names = 2)
pnames <- data.frame(PaperName = character(0), Group = character(0))
yamin <- yaml.load_file(args$yaml)
pheno_names <- c()
for (n in names(yamin$phenotypes)) {
  pheno_names <- c(pheno_names, n)
  pnames <-
    rbind(
      pnames,
      list(
        Group = yamin$phenotypes[[n]]$group,
        PaperName = yamin$phenotypes[[n]]$papername
      )
    )
}
row.names(pnames) <- pheno_names

phenos <- as.character(pnames[phenos, "PaperName", drop = T])
group_phenos <-
  # SNPs annotations
  anno <-
  read_delim(
    paste0(args$outdir, "/annotations.csv"),
    ",",
    col_names = c("rs", "ps", "chr"),
    guess_max = Inf
  )
# Read the p-values from all the pasted phenotypes and mv LOCO files
pvalmat <- NULL

# Read the genotypes
geno_t <-
  read_csv(
    paste0(args$outdir, "/strains_genotypes_all.csv"),
    col_types = cols(
      .default = col_double(),
      chr = col_character(),
      rs = col_character(),
      major = col_character(),
      minor = col_character()
    )
  )

geno <-
  as.matrix(
    geno_t %>% column_to_rownames(var = "rs") %>% dplyr::select(-chr, -bp38, -major, -minor)
  )

PVE <- read_csv(paste0(args$outdir, "/PVE_GEMMA_estimates.txt"))
PVE <-
  left_join(PVE, as_tibble(pnames, rownames = "phenotype"), by = ("phenotype"))
pnames <-
  left_join(pnames, tibble(Group = unique(pnames$Group), color = grpcol[1:length(unique(pnames$Group))]), by =
              "Group")

# We're all set
dir.create(args$plotdir, recursive = TRUE)
set.seed(490)

# Write the genes for INRICH, get the gene annotations and pass them on
annot <- write_genes_map(args$plotdir)

# Plot each phenotype's Manhattan plot
lilp <- vector("list", 0)
allpeaks <- NULL
allsnps <- NULL
all_ispeak <- NULL
all_choose <- NULL
for (i in 1:length(phenos)) {
  pp <-
    plot_gemma_lmm(
      paste0(args$outdir, "/output/lmm_pheno_", i, "_all_LOCO.assoc.txt"),
      name = "Chromosome",
      genotypes = geno,
      namethr = args$pvalthr,
      redthr = args$pvalthr,
      maxdist = args$ldpeakdist,
      corrthr = args$peakcorr,
      annot = annot
    )
  pname <- phenos[i]
  lilp[[pname]] <- pp
  allsnps <- pp$gwas$rs
  if (is.null(pvalmat)) {
    pvalmat <-
      pp$gwas %>% mutate(!!(pname) := P) %>% dplyr::select(rs, !!(pname))
    all_ispeak <-
      pp$gwas %>% mutate(!!(pname) := ispeak) %>% dplyr::select(rs, !!(pname))
    all_choose <-
      pp$gwas %>% mutate(!!(pname) := choose) %>% dplyr::select(rs, !!(pname))
  } else{
    pvalmat <-
      left_join(pvalmat,
                pp$gwas %>% mutate(!!(pname) := P) %>% dplyr::select(rs, !!(pname)),
                by = "rs")
    all_ispeak <-
      left_join(all_ispeak,
                pp$gwas %>% mutate(!!(pname) := ispeak) %>% dplyr::select(rs, !!(pname)),
                by = "rs")
    all_choose <-
      left_join(all_choose,
                pp$gwas %>% mutate(!!(pname) := choose) %>% dplyr::select(rs, !!(pname)),
                by = "rs")
  }
  allpeaks <- c(allpeaks, pp$gwas$rs[pp$gwas$ispeak])
  ggsave(
    filename = paste0(
      args$plotdir,
      "/Manhattan_plot_phenotype_",
      i,
      "_",
      phenos[i],
      ".pdf"
    ),
    plot = pp$plot + theme(text = element_text(size = 10, family =
                                                 ffam)),
    dpi = "print",
    device = cairo_pdf,
    width = fullw,
    height = height,
    units = "in"
  )
}
# Plot the groups MAnhattan plots
grpwas <- list()
if (args$nomv) {
  # Plot each group's max P
  for (g in c(unique(as.character(pnames$Group)), "All Phenotypes")) {
    allpwas = NULL
    plist <- pnames$PaperName[pnames$Group == g]
    if (g == "All Phenotypes") {
      plist <- pnames$PaperName
    }
    for (p in intersect(plist, names(lilp))) {
      if (is.null(allpwas)) {
        allpwas <- lilp[[p]]$pwas %>% dplyr::select(-ispeak,-choose)
      } else{
        allpwas <-
          left_join(allpwas,
                    lilp[[p]]$pwas[, c("rs", "P", "p_wald")],
                    by = "rs",
                    suffix = c("", ".x"))
        allpwas$P <- pmax(allpwas$P, allpwas$P.x)
        allpwas$p_wald <- pmin(allpwas$p_wald, allpwas$p_wald.x)
        allpwas <- allpwas %>% dplyr::select(-P.x,-p_wald.x)
      }
    }
    if (is.null(allpwas)) {
      print(g)
      next
    }
    pnums <-
      rep_peaks(
        geno,
        allpwas,
        rs_thr = args$peakcorr,
        pthr = 10 ^ -args$pvalthr,
        mxd = args$ldpeakdist
      )
    allpwas <- allpwas %>% left_join(pnums, by = "rs")
    pname = g
    allpeaks <- c(allpeaks, allpwas$rs[allpwas$ispeak])
#   pvalmat <-
#      left_join(pvalmat,
#                allpwas %>% mutate(!!(pname) := P) %>% dplyr::select(rs, !!(pname)),
#                by = "rs")

    grpwas[[g]] <- allpwas
    # Recolor the second layer with the clusters colors
  }
} else{
  for (grpf in Sys.glob(paste0(
    args$outdir,
    "/output/lmm_phenotypes_*_all_LOCO.assoc.txt"
  ))) {
    pp <-
      plot_gemma_lmm(
        grpf,
        name = "Chromosome",
        genotypes = geno,
        namethr = args$pvalthr,
        redthr = args$pvalthr,
        maxdist = args$lgpeakdist,
        corrthr = args$peakcorr,
        test = "p_score",
        annot = annot
      )
    if (is.null(mp))
      mp <- pp$plot
    pname <-
      gsub(".*phenotypes_(.*)_all_LOCO.assoc.txt", "\\1", grpf)
    grpwas[[pname]] = pp$pwas

    allpeaks <- c(allpeaks, pp$gwas$rs[pp$gwas$ispeak])
    pvalmat <-
      left_join(pvalmat,
                pp$gwas %>% mutate(!!(pname) := P) %>% dplyr::select(rs, !!(pname)),
                by = "rs")

    ggsave(
      filename = paste0(
        args$plotdir,
        "/Manhattan_plot_phenotypes_",
        pname,
        ".pdf"
      ),
      plot = pp$plot + theme(text = element_text(size = 10, family =
                                                   ffam)),
      dpi = "print",
      device = cairo_pdf,
      width = fullw,
      height = height,
      units = "in"
    )
  }
}

# Cluster the peaks using the P values

pgwas <-
  pvalmat %>% filter(rs %in% allpeaks) %>% column_to_rownames(var = "rs")
print(sum(is.na(pgwas)))
pgwas[is.na(pgwas)] = 0
pgwas <- as.matrix(pgwas)

kk <- kmeans(pgwas, args$clusters)

# Plot the m-value heatmap
# Order the columns
rowarr <- NULL
for (i in 1:args$cluster) {
  cc <- hclust(dist(pgwas[kk$cluster == i, , drop = F]))
  rowarr <- c(rowarr, row.names(pgwas)[kk$cluster == i][cc$order])
}
pgwas <- pgwas[rowarr,]
clustcol <- tibble(cluster = 1:args$clusters, color = ccols)
colrow <-
  tibble(rs = rownames(pgwas), cluster = kk$cluster[rowarr]) %>% left_join(clustcol, by =
                                                                             "cluster") %>% column_to_rownames(var = "rs") %>% dplyr::select(color)
grptocol <-
  tibble(Group = unique(PVE$Group), color = grpcol[1:length(unique(PVE$Group))])
colcol <-
  PVE %>% left_join(grptocol) %>% column_to_rownames(var = "PaperName") %>% dplyr::select(color)
colcol <- colcol$color

cairo_pdf(
  paste0(args$plotdir, "/all_peaks_heatmap.pdf"),
  width = fullw,
  height = height + 1,
  family = ffam
)
heatmap.2(
  pgwas,
  col = hmcol,
  Rowv = F,
  Colv = T,
  dendrogram = "col",
  scale = "none",
  trace = "none",
  RowSideColors = colrow[, 1, drop = T],
  ColSideColors = colcol,
  labRow = NA,
  hclustfun = function(x)
    hclust(x, method = "average"),
  distfun = function(x)
    dist(scale(x)),
  margins = c(12, 8),
  srtCol = 45,
  key = T,
  density.info = "none"
)
dev.off()


# Plot the PVE estimates with SE
pvep <-
  ggplot(PVE, aes(reorder(PaperName, -PVE), PVE, fill = Group)) + geom_bar(color =
                                                                             "black", stat = "identity") +
  scale_fill_manual(values = grpcol) +
  geom_errorbar(aes(ymin = PVE - PVESE, ymax = PVE + PVESE), width = .2) +
  xlab("Phenotype") +
  theme_bw() + theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  )) +
  theme(text = element_text(size = 10, family = ffam))
ggsave(
  paste0(args$plotdir, "/PVE_plot.pdf"),
  plot = pvep,
  device = cairo_pdf,
  dpi = "print",
  width = fullw,
  height = height,
  units = "in"
)


# Plot the phenotypes manhattan plots with clusters colors
# Add the cluster number to the pwas object
for (i in names(lilp)) {
  p <- lilp[[i]]
  p$pwas <-
    p$pwas %>% left_join(tibble(rs = rownames(pgwas), cluster = as.factor(kk$cluster[rowarr])), by =
                           "rs")
  if (sum(p$pwas$ispeak) == 0)
    next
  p$pwas <-
    p$pwas %>% dplyr::select(-cluster) %>% left_join(filter(p$pwas, ispeak) %>%
                                                       dplyr::select(choose, cluster),
                                                     by = "choose")

  # Recolor the second layer with the clusters colors
  pnoname <- p$plot
  pnoname$layers <- pnoname$layers[1:2]
  ggsave(
    filename = paste0(args$plotdir, "/replot_Manhattan_clusters_", i, ".pdf"),
    plot = pnoname + ggnewscale::new_scale_color() +
      geom_point(aes(color = p$pwas$cluster), size = 0.9) +
      scale_color_manual(values = ccols) +
      ggnewscale::new_scale_color() +
      geom_point(
        aes(alpha = p$pwas$ispeak),
        size = 1.2,
        color = "black"
      ) +
      scale_alpha_manual(values = c(0, 1)) +
      ggnewscale::new_scale_color() +
      geom_point(aes(
        color = p$pwas$cluster,
        alpha = p$pwas$ispeak
      ), size = 0.9) +
      scale_color_manual(values = ccols) +
      scale_alpha_manual(values = c(0, 1)) +
      theme(text = element_text(size = 10, family = ffam)),
    device = cairo_pdf,
    dpi = "print",
    width = fullw,
    height = height,
    units = "in"
  )
}


for (g in names(grpwas)) {
  # Add cluster to ispeak
  allpwas <-
    grpwas[[g]] %>% left_join(tibble(rs = rownames(pgwas), cluster = as.factor(kk$cluster[rowarr])), by =
                                "rs")
  # Expand cluster to all choose
  allpwas <-
    allpwas %>% dplyr::select(-cluster) %>% left_join(filter(allpwas, ispeak) %>%
                                                        dplyr::select(choose, cluster),
                                                      by = "choose")
  axisdf <-
    allpwas %>% group_by(chr) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ymax <- 1.25 * max(allpwas$P, na.rm = TRUE)
  ymin <- 1.25 * min(allpwas$P, na.rm = TRUE)
  chr_label <- axisdf$chr
  chr_label[chr_label == 20] = "X"
  ggsave(
    filename = paste0(
      args$plotdir,
      "/replot_Manhattan_clusters_",
      gsub(" ", "_", g),
      ".pdf"
    ),
    plot = ggplot2::ggplot(allpwas, aes(x = BPcum, y = P)) +

      # Show all points
      geom_point(aes(color = as.factor(chr)) , alpha = 1, size = 0.7) +
      scale_color_manual(values = c(rep(
        c("#CCCCCC", "#969696"), 10
      ))) +
      ggnewscale::new_scale_color() +
      geom_point(aes(color = cluster), size = 0.9) +
      scale_color_manual(values = ccols) +
      ggnewscale::new_scale_color() +
      geom_point(aes(alpha = ispeak), size = 1.2, color = "black") +
      scale_alpha_manual(values = c(0, 1)) +
      ggnewscale::new_scale_color() +
      geom_point(aes(color = cluster,
                     alpha = ispeak), size = 0.9) +
      scale_color_manual(values = ccols) +
      scale_alpha_manual(values = c(0, 1)) +
      scale_x_continuous(label = chr_label, breaks = axisdf$center) +
      scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
      ylim(ymin, ymax) +
      xlab(g) +
      ylab("-log(P-value)") +
      theme_bw() +
      theme(
        legend.position = "none",
        text = element_text(size = 10, family = ffam),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      ),
    device = cairo_pdf,
    dpi = "print",
    width = fullw,
    height = height,
    units = "in"
  )
}


# Plot the LD drop figure
comp_LD_2 <-
  function(genotypes,
           c,
           maxdist = 2500000,
           MAF = 0.1,
           miss = 0.1) {
    gen <- genotypes %>% filter(chr == c, !is.na(rs)) %>% arrange(bp38)
    gen <-
      gen[(
        rowSums(gen[, 6:ncol(gen)] == 0, na.rm = T) > (dim(gen)[2] - 5) * MAF &
          rowSums(gen[, 6:ncol(gen)] > 0, na.rm = T) > (dim(gen)[2] -
                                                          5) * MAF &
          rowSums(is.na(gen[, 6:ncol(gen)])) < (dim(gen)[2] - 5) *
          miss
      ),]
    allcor = NULL
    for (i in 1:(nrow(gen) - 1)) {
      tmat <-
        as.data.frame(base::t(
          gen %>% filter(between(bp38, gen$bp38[i], gen$bp38[i] + maxdist)) %>%
            dplyr::select(-bp38, -chr, -major, -minor) %>% column_to_rownames(var = "rs")
        ))
      if (ncol(tmat) < 2)
        next()
      ct <-
        cor(tmat[, 1], tmat[, 2:ncol(tmat), drop = F], use = 'pairwise.complete.obs', method =
              "pearson")
      allcor <-
        rbind(
          allcor,
          as_tibble(t(ct), rownames = "SNP2") %>% mutate(r_sq = V1 ^ 2, SNP1 = gen$rs[i]) %>% dplyr::select(SNP1, SNP2, r_sq)
        )
    }
    return(allcor)
  }
geno_s <-
  geno_t %>% filter(chr != "X", chr != "Y", chr != "MT") %>% sample_n(args$sample)
allchr = NULL
for (chr in unique(geno_s$chr)) {
  allchr <- rbind(allchr, comp_LD_2(geno_s, chr))
}
allchr <-
  allchr %>% left_join(dplyr::select(geno_s, rs, bp38), by = c("SNP1" = "rs")) %>%
  left_join(dplyr::select(geno_s, rs, bp38), by = c("SNP2" = "rs")) %>% mutate(dist = bp38.y -
                                                                                 bp38.x)
avgwin = 5000
allchr$distc <-
  (cut(allchr$dist, breaks = seq(
    from = min(allchr$dist) - 1,
    to = max(allchr$dist) + 1,
    by = avgwin
  )))
allavg <-
  allchr %>% group_by(distc) %>% summarise(avdist = mean(dist), avr_sq =
                                             mean(r_sq, na.rm = T)) %>% ungroup()
pld <-
  ggplot(allavg, aes(avdist / 1000000, avr_sq)) + geom_smooth(
    color = RColorBrewer::brewer.pal(3, "Set1")[3],
    method = "loess",
    span = 0.3,
    se = FALSE
  ) +
  xlim(c(0, 2.5)) + labs(x = "Distance (Mbp)", y = expression("Average LD" ~
                                                                (r ^ {
                                                                  2
                                                                })))
ggsave(
  filename = paste0(args$plotdir, "/plot_LD_drop.pdf"),
  plot = pld + theme_bw() + theme(
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(size = 10, family = ffam)
  ),
  device = cairo_pdf,
  dpi = "print",
  width = halfw,
  height = height,
  units = "in"
)


# Plot MAF histogram
mafdat <-
  tibble(rs = geno_t$rs, maf = rowSums(geno_t[, -1:-5]) / (2 * (ncol(geno_t) -
                                                                  5)))
mafdat$maf <- pmin(mafdat$maf, 1 - mafdat$maf)
mafdat <- mafdat %>% filter(rs %in% allsnps)
if ("All Phenotypes" %in% names(grpwas)) {
  mafdat <- inner_join(mafdat, grpwas[["All Phenotypes"]] %>% dplyr::select(rs, choose), by = "rs")
} else{
  mafdat$choose <- 0
}
mafp <-
  ggplot(mafdat, aes(maf, fill = (choose > 0), color = (choose > 0))) + geom_histogram(binwidth = 1 /
                                                                                       (ncol(geno_t) - 5)) + xlim(c(0, 0.5)) +
  scale_color_manual(
    values = RColorBrewer::brewer.pal(12, "Paired")[3:4],
    name = "",
    labels = c("All", "Peak")
  ) +
  scale_fill_manual(
    values = RColorBrewer::brewer.pal(12, "Paired")[3:4],
    name = "",
    labels = c("All", "Peak")
  ) +
  xlab("MAF") +
  theme_bw() + theme(
    legend.position = c(0.8, 0.9),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(size = 10, family = ffam)
  )
ggsave(
  filename = paste0(args$plotdir, "/plot_MAF_hist.pdf"),
  plot = mafp,
  device = cairo_pdf,
  dpi = "print",
  width = halfw,
  height = height,
  units = "in"
)




# Run enrichR for each phenotype
dbs <- listEnrichrDbs()

# Expand each peak to include the entire peak, not just the single SNP
ext_peak_sing <- function(snps, maxdist = 500000) {
  csum <-
    snps %>% group_by(choose) %>% dplyr::summarize(maxps = max(ps), minps = min(ps)) %>% ungroup()
  snps %>% left_join(csum, by = "choose") %>% mutate(maxps = pmin(maxps, ps +
                                                                    maxdist),
                                                     minps = pmax(minps, ps - maxdist))
  #snps %>% mutate(minps = ps - maxdist, maxps = ps + maxdist)
}

# This tibble will accumulate all the genes for each cluster (ID and name)
clustgene <- vector(mode = "list", length = args$clusters)
clustmgi <- vector(mode = "list", length = args$clusters)
allgenes = c()
pgc <- tibble(rs = rownames(pgwas), cluster = kk$cluster[rowarr])
clusterpeaks = tibble(
  cluster = numeric(0),
  chr = character(0),
  minps = numeric(0),
  maxps = numeric(0)
)
write_inrich_snps(geno_t, args$plotdir)
for (i in 1:length(lilp)) {
  pp <- lilp[[i]]
  if (sum(pp$gwas$ispeak) == 0)
    next
  expp <- ext_peak_sing(pp$gwas, maxdist = args$peakdist)
  write_inrich_phenotype(expp[expp$ispeak == T, ], args$plotdir, phenos[i])
  run_inrich(
    args$plotdir,
    phenos[i],
    exec = args$inrich,
    i = args$inrich_i,
    j = args$inrich_j
  )
  exppc <-
    expp %>% dplyr::filter(ispeak) %>% left_join(pgc, by = "rs") %>% dplyr::select(cluster, chr, minps, maxps)
  clusterpeaks <- rbind(clusterpeaks, exppc)
  affgen <-
    get_genes(expp[expp$ispeak == T, ], dist = 1000, annot = annot)
  if (nrow(affgen) > 0) {
    # Add the genes to the appropriate cluster
    clustp <- left_join(affgen, pgc, by = "rs")
    for (c in unique(clustp$cluster)) {
      clustgene[[c]] <-
        unique(c(clustgene[[c]], clustp$ensembl_gene_id[clustp$cluster == c]))
      clustmgi[[c]] <-
        unique(c(clustmgi[[c]], clustp$mgi_symbol[clustp$cluster == c]))
    }
    allgenes <- unique(c(allgenes, clustp$mgi_symbol))
    write_csv(affgen,
              path = paste0(args$plotdir, "/genes_for_phenotype_", i, ".csv"))
    # Run enrichr
    enrr <-
      enrichr(unique(affgen$mgi_symbol[!(
        grepl(pattern = "^Gm", x =  affgen$mgi_symbol) |
          grepl("Rik$", affgen$mgi_symbol) |
          affgen$mgi_symbol == ""
      )]), dbs$libraryName)
    for (d in dbs$libraryName) {
      if (length(dim(enrr[[d]])) > 1 && dim(enrr[[d]])[2] > 1) {
        rtb <- as_tibble(enrr[[d]]) %>% filter(Adjusted.P.value <= 0.05)
        if (nrow(rtb) > 0)
          write_csv((
            rtb %>% mutate(
              total.genes = length(affgen$mgi_symbol[!(
                grepl(pattern = "^Gm", x =  affgen$mgi_symbol) |
                  grepl("Rik$", affgen$mgi_symbol) |
                  affgen$mgi_symbol == ""
              )]),
              phenotype = pnames$PaperName[i],
              library = d
            )
          ),
          path = paste0(
            args$plotdir,
            "/enrichR_phenotype_",
            i,
            "_db_",
            d,
            "p005.csv"
          )
          )

      }
    }
  }
}

# Run each group
for (n in names(grpwas)) {
  pp <- grpwas[[n]]
  # Change the chromosome to character
  pp$chr <- as.character(pp$chr)
  pp$chr[pp$chr == "20"] <- "X"
  if (sum(pp$ispeak) == 0)
    next
  expp <- ext_peak_sing(pp, maxdist = args$peakdist)
  write_inrich_phenotype(expp[expp$ispeak == T, ], args$plotdir, n)
  run_inrich(
    args$plotdir,
    n,
    exec = args$inrich,
    i = args$inrich_i,
    j = args$inrich_j
  )
  affgen <-
    get_genes(expp[expp$ispeak == T, ], dist = 1000, annot = annot)
  if (nrow(affgen) > 0) {
    # Add the genes to the appropriate cluster
    write_csv(affgen,
              path = paste0(
                args$plotdir,
                "/genes_for_phenotype_Group_",
                gsub(" ", "_", n),
                ".csv"
              ))
    # Run enrichr
    enrr <-
      enrichr(unique(affgen$mgi_symbol[!(
        grepl(pattern = "^Gm", x =  affgen$mgi_symbol) |
          grepl("Rik$", affgen$mgi_symbol) |
          affgen$mgi_symbol == ""
      )]), dbs$libraryName)
    for (d in dbs$libraryName) {
      if (length(dim(enrr[[d]])) > 1 && dim(enrr[[d]])[2] > 1) {
        rtb <- as_tibble(enrr[[d]]) %>% filter(Adjusted.P.value <= 0.05)
        if (nrow(rtb) > 0)
          write_csv((
            rtb %>% mutate(
              total.genes = length(affgen$mgi_symbol[!(
                grepl(pattern = "^Gm", x =  affgen$mgi_symbol) |
                  grepl("Rik$", affgen$mgi_symbol) |
                  affgen$mgi_symbol == ""
              )]),
              phenotype = n,
              library = d
            )
          ),
          path = paste0(
            args$plotdir,
            "/enrichR_phenotype_Group_",
            gsub(" ", "_", n),
            "_db_",
            d,
            "p005.csv"
          )
          )

      }
    }
  }
}

for (k in 1:args$clusters) {
  write_csv(
    data.frame(genes = clustgene[[k]]),
    path = paste0(args$plotdir, "/genes_for_cluster_", k, ".csv")
  )
  # Run INRICH
  write_inrich_phenotype(clusterpeaks %>% filter(cluster == k),
                         args$plotdir,
                         paste0("cluster_", k))
  run_inrich(
    args$plotdir,
    paste0("cluster_", k),
    exec = args$inrich,
    i = args$inrich_i,
    j = args$inrich_j
  )

  # Run enrichr
  glen = length(clustmgi[[k]][!(grepl(pattern = "^Gm", x = clustmgi[[k]]) |
                                  grepl("Rik$", clustmgi[[k]]))])
  enrr <-
    enrichr(unique(clustmgi[[k]][!(grepl(pattern = "^Gm", x = clustmgi[[k]]) |
                                     grepl("Rik$", clustmgi[[k]]))]), dbs$libraryName)
  for (d in dbs$libraryName) {
    if (dim(enrr[[d]])[2] > 1) {
      rtb <- as_tibble(enrr[[d]]) %>% filter(Adjusted.P.value <= 0.05)
      if (nrow(rtb) > 0)
        write_csv((rtb %>% mutate(
          total.genes = glen,
          cluster = k,
          library = d
        )),
        path = paste0(
          args$plotdir,
          "/enrichR_cluster_",
          k,
          "_db_",
          d,
          "p005.csv"
        )
        )

    }
  }
}

# Run INRICH
write_inrich_phenotype(clusterpeaks, args$plotdir, "all_phenotypes")
run_inrich(
  args$plotdir,
  "all_phenotypes",
  exec = args$inrich,
  i = args$inrich_i,
  j = args$inrich_j
)

# Plot markers density
chrord <- c("X", 19:1)
densp <- geno_t %>% filter(chr != "Y", chr != "MT")  %>%
  ggplot(aes(bp38 / 1000000, factor(chr, levels = chrord))) +
  geom_bin2d(binwidth = 1, drop = T) + xlab("Position (Mbp)") + ylab ("Chromosome") +
  scale_fill_viridis(name = expression(frac('markers', '1 Mbp'))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(size = 10, family = ffam)
  )
ggsave(
  filename = paste0(args$plotdir, "/Chromosome_density_plot.pdf"),
  plot = densp,
  device = cairo_pdf,
  dpi = "print",
  width = fullw,
  height = height,
  units = "in"
)
