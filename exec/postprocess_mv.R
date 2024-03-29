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
library(yaml)
library(cowplot)
library(grid)
library(gtable)
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
parser$add_argument("--ldpeakdist",
                    type = "integer",
                    default = 10000000,
                    help = "Peak maximal width")
parser$add_argument("--peakcorr",
                    type = "double",
                    default = 0.25,
                    help = "SNPs R^2 correlation cutoff for peak determination")
parser$add_argument("--loddrop",
                    type = "double",
                    default = 1.5,
                    help = "LOD drop from the peak SNP for the purpose of selecting genes related to the QTL")
parser$add_argument(
  "--external_inrich",
  action = "store_true",
  default = FALSE,
  help = "Prepare INRICH files but don't run INRICH"
)
parser$add_argument(
  "--coat_phenotype",
  action = "store_true",
  default = FALSE,
  help = "GWAS of coat color, no phenotypes in yaml file"
)
parser$add_argument("--colorgroup",
                    action = "store_true",
                    default = FALSE,
                    help = "Color the Group Manhattan plots by group rather than cluster")
parser$add_argument("--meanvariance", action="store_true", default=FALSE,
                    help="Plot the PVE of mean and variance phenotypes in different plots")
parser$add_argument("--set3", action="store_true", default=FALSE,
                    help="Use Set3 color palette for groups, default is Accent")
parser$add_argument("--minherit", type="double", default=0,
                    help= "Heritability threshold (PVE) to include a phenotype in the aggregated reports")
parser$add_argument("--ploteffect", action="store_true", default=FALSE, help="Plot effect plots for all QTL peaks")
args <- parser$parse_args()

# Step 1: Read the color palette
# Cluster colors
ccols <- brewer.pal(args$clusters, "Dark2")[1:args$clusters]
# Heatmap plot for m-values
pigr <- RColorBrewer::brewer.pal(name = "PiYG", n = 11)
hmcol <-
  viridis(128)#colorRampPalette(pigr[c(2,5,10)])(128)#viridis(128, option="cividis")
if (args$set3) {
  grpcol <- RColorBrewer::brewer.pal(12, "Set3")[2:12]
} else{
  grpcol <- RColorBrewer::brewer.pal(8, "Accent")
}
fullw <- 7.25
halfw <- 3.54
fheight <- 11-1.25
height <- 3.54#fheight/4

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
pnames <- data.frame(PaperName = character(0), Group = character(0), stringsAsFactors = FALSE)
yamin <- yaml.load_file(args$yaml)
pheno_names <- c()
for (n in names(yamin$phenotypes)) {
  pheno_names <- c(pheno_names, n)
  pnames <-
    rbind(
      pnames,
      data.frame(
        Group = if (length(yamin$phenotypes[[n]]$group)) yamin$phenotypes[[n]]$group else "NoGroup",
        PaperName = yamin$phenotypes[[n]]$papername, stringsAsFactors = FALSE
      )
    )
}
if (args$coat_phenotype){
  for (p in phenos){
    pheno_names <- c(pheno_names, p)
    pnames <- rbind(pnames, data.frame(Group=if(grepl("coat", p)) "Coat" else "Eyes", PaperName=gsub(" ", "-", gsub("(coat|eyes)", "", p))))
  }
}
row.names(pnames) <- pheno_names

phenos <- as.character(pnames[phenos, "PaperName", drop = T])

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
groupsOrder = if (length(yamin$groups)) yamin$groups else unique(pnames$Group)
PVE <-
  left_join(PVE, as_tibble(pnames, rownames = "phenotype"), by = ("phenotype"))
if (args$meanvariance){
  gnames <- groupsOrder[!grepl("Variance", groupsOrder)]
  grpcol <- rep(grpcol[1:min(length(gnames), length(grpcol))], ceiling(length(gnames)/length(grpcol)))
  grpcol <- grpcol[1:length(gnames)]
  names(grpcol) <- gnames
  for (g in groupsOrder[grepl("Variance", groupsOrder)]){
    cadd <- grpcol[gsub("Variance", "", g)]
    names(cadd) <- g
    grpcol <- c(grpcol, cadd)
  }
}else{
grpcol <- rep(grpcol[1:min(length(groupsOrder), length(grpcol))], ceiling(length(groupsOrder)/length(grpcol)))
grpcol <- grpcol[1:length(groupsOrder)]
names(grpcol) <- groupsOrder
}
pnames <-
  left_join(pnames, tibble(Group = groupsOrder, color = grpcol), by =
              "Group")
row.names(pnames) <- pheno_names

# We're all set
dir.create(args$plotdir, recursive = TRUE)
set.seed(490)

# Plot the PVE estimates with SE

pvh <- height
if (args$meanvariance) {
  if (dim(PVE)[1] > 40)
    pvh <- height * 2
  pveplot <- paired_PVE_plot(PVE, minherit=args$minherit)
  pvep <-
    cbind(
      ggplotGrob(pveplot$mean_plot + scale_fill_manual(values = grpcol) + theme_bw()  +
        theme(
          text = element_text(size = 10, family = ffam),
          legend.position = "none",
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )),
      ggplotGrob(pveplot$var_plot + scale_fill_manual(values = grpcol) + theme_bw()  +
        theme(
          text = element_text(size = 10, family = ffam),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
        ))
    )
  cairo_pdf(
    paste0(args$plotdir, "/PVE_plot.pdf"),
    width = halfw * 1.5,
    height = pvh,
    family = ffam
  )
  grid.draw(pvep)
  dev.off()
} else{
  if (dim(PVE)[1] > 40)
    pvh <- height * 2
  pvep <-
    ggplot(PVE, aes(reorder(PaperName,-PVE), PVE, fill = Group, alpha=PVE>args$minherit)) + geom_bar(color =
                                                                              "black", stat = "identity") +
    scale_fill_manual(values = grpcol) +
    geom_errorbar(aes(ymin = PVE - PVESE, ymax = PVE + PVESE), width = .2) +
    scale_alpha_manual(breaks = c(FALSE, TRUE), values=c(0.4,1), guide=FALSE) +
    xlab("Phenotype") + coord_flip() +
    theme_bw()  +
    theme(text = element_text(size = 10, family = ffam),
          legend.position = "right")

ggsave(
  paste0(args$plotdir, "/PVE_plot.pdf"),
  plot = pvep,
  device = cairo_pdf,
  dpi = "print",
  width = halfw * 1.5,
  height = pvh,
  units = "in"
)
}

# Write the genes for INRICH, get the gene annotations and pass them on
annot <- write_genes_map(args$plotdir)

# Plot each phenotype's Manhattan plot
lilp <- vector("list", 0)
allpeaks <- NULL
allsnps <- NULL

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
  } else{
    pvalmat <-
      left_join(pvalmat,
                pp$gwas %>% mutate(!!(pname) := P) %>% dplyr::select(rs, !!(pname)),
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
    plot = pp$plot + ggtitle(pname) + theme(text = element_text(size = 10, family =
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
  # Add "All Phenotypes" to the list of groups. If meanvariance add two groups for the mean and the variance
  grplist <- c(unique(as.character(pnames$Group)), "All Phenotypes")
  if (args$meanvariance){
    grplist <- c(grplist, "All Mean Phenotypes", "All Variance Phenotypes")
  }
  for (g in grplist) {
    allpwas = NULL
    plist <- pnames$PaperName[pnames$Group == g]
    if (g == "All Phenotypes") {
      plist <- pnames$PaperName
    }
    if (g == "All Mean Phenotypes"){
      plist <- pnames$PaperName[!grepl("Variance", pnames$PaperName)]
    }
    if (g == "All Variance Phenotypes"){
      plist <- pnames$PaperName[grepl("Variance", pnames$PaperName)]
    }
    plist <- intersect(plist, PVE$PaperName[PVE$PVE >= args$minherit])
    for (p in intersect(plist, names(lilp))) {
      if (is.null(allpwas)) {
        allpwas <- lilp[[p]]$pwas %>% dplyr::select(chr, rs, ps, allele0, allele1, af, P, p_wald, BPcum)
        allpwas$grpcolor <- pnames$Group[pnames$PaperName==p][1]
      } else{
        allpwas <-
          left_join(allpwas,
                    lilp[[p]]$pwas[, c("rs", "P", "p_wald")],
                    by = "rs",
                    suffix = c("", ".x"))
        allpwas$P <- pmax(allpwas$P, allpwas$P.x)
        allpwas$p_wald <- pmin(allpwas$p_wald, allpwas$p_wald.x)
        allpwas$grpcolor[allpwas$p_wald == allpwas$p_wald.x] <- pnames$Group[pnames$PaperName==p][1]
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
        maxdist = args$ldpeakdist,
        corrthr = args$peakcorr,
        test = "p_score",
        annot = annot
      )
    pname <-
      gsub(".*phenotypes_(.*)_all_LOCO.assoc.txt", "\\1", grpf)
    grpwas[[pname]] = pp$pwas

    allpeaks <- c(allpeaks, pp$gwas$rs[pp$gwas$ispeak])
#    pvalmat <-
#      left_join(pvalmat,
#                pp$gwas %>% mutate(!!(pname) := P) %>% dplyr::select(rs, !!(pname)),
#                by = "rs")

    ggsave(
      filename = paste0(
        args$plotdir,
        "/Manhattan_plot_phenotypes_",
        pname,
        ".pdf"
      ),
      plot = pp$plot + ggtitle(pname) + theme(text = element_text(size = 10, family =
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
for (i in unique(kk$cluster)) {
  cc <- hclust(dist(pgwas[kk$cluster == i, , drop = F]))
  rowarr <- c(rowarr, row.names(pgwas)[kk$cluster == i][cc$order])
}
pgwas <- pgwas[rowarr,]
clustcol <- tibble(cluster = unique(kk$cluster), color = ccols)
colrow <-
  tibble(rs = rownames(pgwas), cluster = kk$cluster[rowarr]) %>% left_join(clustcol, by =
                                                                             "cluster") %>% column_to_rownames(var = "rs") %>% dplyr::select(color)
grptocol <-
  tibble(Group = groupsOrder, color = grpcol[1:length(groupsOrder)])
colcol <-
  PVE %>% left_join(grptocol) %>% column_to_rownames(var = "PaperName") %>% dplyr::select(color)
colcol <- colcol$color
hwid <- fullw
if (dim(PVE)[1]>40) hwid <- fullw*2

write.csv(cbind(pgwas, cbind(colrow, kk$cluster[rowarr])), file = paste0(args$plotdir, "/all_peaks_values_clusters.csv"))
cairo_pdf(
  paste0(args$plotdir, "/all_peaks_heatmap.pdf"),
  width = hwid,
  height = height + 1,
  family = ffam
)
hplt <- heatmap.2(
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
  key.title = "-log(P-value)",
  density.info = "none"
)

dev.off()

svg(paste0(args$plotdir, "/all_peaks_heatmap.svg"),
    width = hwid,
    height = height + 1,
    family = ffam)
eval(hplt$call)
dev.off()


# Plot the phenotypes manhattan plots with clusters colors
# Add the cluster number to the pwas object
if (!args$colorgroup){
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
    plot = pnoname + ggnewscale::new_scale("alpha") + ggnewscale::new_scale("color") + ggnewscale::new_scale("size") +
      geom_point(aes(color = p$pwas$cluster, size=P, alpha=rsq)) +
      scale_color_manual(values = ccols) +
      scale_size_continuous(range=c(0,1), trans = "exp") +
      scale_alpha_continuous(range = c(0,1), trans="exp") +
      ggnewscale::new_scale("alpha") + ggnewscale::new_scale("color") + ggnewscale::new_scale("size") +
      geom_point(
        aes(alpha = p$pwas$ispeak),
        size = 1.2,
        color = "black"
      ) +
      scale_alpha_manual(values = c(0, 1)) + ggnewscale::new_scale("alpha") + ggnewscale::new_scale("color") + ggnewscale::new_scale("size") +
      #ggnewscale::new_scale_color() +
      geom_point(aes(
        color = p$pwas$cluster,
        alpha = p$pwas$ispeak
      ), size = 1) +
      scale_color_manual(values = ccols) +
      scale_alpha_manual(values = c(0, 1)) +
      ggtitle(i) +
      theme(text = element_text(size = 10, family = ffam)),
    device = cairo_pdf,
    dpi = "print",
    width = fullw,
    height = height,
    units = "in"
  )
}
}
mainplot = NULL
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
  colorby <- "cluster"
  palette <- ccols
  pbreaks <- unique(kk$cluster)
  if (args$colorgroup) {
    colorby = "grpcolor"
    palette <- grpcol
  }
  print(palette)
  outplot <- ggplot2::ggplot(allpwas, aes(x = BPcum, y = P)) +

    # Show all points
    geom_point(aes(color = as.factor(chr), size=P) , alpha = 1) +
    scale_color_manual(values = c(rep(
      c("#CCCCCC", "#969696"), 10
    ))) +
    scale_size_continuous(range=c(0,1), trans = "exp") +
    geom_segment(y = args$pvalthr, x=min(allpwas$BPcum)-50000000, xend=max(allpwas$BPcum)+50000000, yend=args$pvalthr,color="#FCBBA1") +
    ggnewscale::new_scale("alpha") + ggnewscale::new_scale("color") + ggnewscale::new_scale("size")  +
    geom_point(aes_string(color = colorby, size="rsq", alpha="rsq")) +
    scale_color_manual(values = palette) +
    scale_size_continuous(range=c(0,1), trans = "exp") +
    scale_alpha_continuous(range = c(0,1), trans="exp") +
    ggnewscale::new_scale("alpha") + ggnewscale::new_scale("color") + ggnewscale::new_scale("size")  +
    geom_point(aes(alpha = ispeak), size = 1.2, color = "black") +
    scale_alpha_manual(values = c(0, 1)) +
    ggnewscale::new_scale("alpha") + ggnewscale::new_scale("color") + ggnewscale::new_scale("size")  +
    geom_point(aes_string(color = colorby,
                   alpha = "ispeak"), size = 1) +
    scale_color_manual(values = palette) +
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
    )
  ggsave(
    filename = paste0(
      args$plotdir,
      "/replot_Manhattan_clusters_",
      gsub(" ", "_", g),
      ".pdf"
    ),
    plot = outplot,
    device = cairo_pdf,
    dpi = "print",
    width = fullw,
    height = height,
    units = "in"
  )
  if (g == "All Phenotypes") mainplot = outplot
  if (g == "All Mean Phenotypes") meanplot = outplot
  if (g == "All Variance Phenotypes") varplot = outplot
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
                                                                }))) + theme_bw() + theme(
                                                                  panel.border = element_blank(),
                                                                  panel.grid.major.x = element_blank(),
                                                                  panel.grid.minor.x = element_blank(),
                                                                  panel.grid.major.y = element_blank(),
                                                                  panel.grid.minor.y = element_blank(),
                                                                  text = element_text(size = 10, family = ffam)
                                                                )
ggsave(
  filename = paste0(args$plotdir, "/plot_LD_drop.pdf"),
  plot = pld,
  device = cairo_pdf,
  dpi = "print",
  width = halfw,
  height = height,
  units = "in"
)
# Plot Figure 1: pvep pld and mainplot
combp <- plot_grid(plot_grid(pvep, pld, ncol=2, nrow=1, labels=c('A', 'B'), label_size = 12, label_fontface = "plain", rel_widths = c(1.5,1)), mainplot, NULL, nrow = 3, ncol = 1, labels = c('', 'C','D'), label_size = 12, fontface="plain")
ggsave(filename = paste0(args$plotdir, "/combined_figure1.svg"),
       plot = combp,
       device = svg,
       dpi = "print",
       width = fullw,
       height = fheight,
       units = "in")
if (args$meanvariance){
  combp <- plot_grid(pvep, meanplot, varplot, mainplot, nrow = 4, ncol = 1, labels = c('A', 'B', 'C', 'D'), label_size = 12, fontface="plain", rel_heights = c(2.5, 1, 1, 1))
  ggsave(filename = paste0(args$plotdir, "/combined_figure2.svg"),
         plot = combp,
         device = svg,
         dpi = "print",
         width = fullw,
         height = fheight,
         units = "in")
}
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


# Expand each peak to include the entire peak, not just the single SNP
ext_peak_sing <- function(snps, maxdist = 500000, loddrop = args$loddrop) {
  if (loddrop){
  snps$minps <- snps$ps
  snps$maxps <- snps$ps
  for (c in unique(snps$choose)){
    peakps <- snps$ps[snps$ispeak & snps$choose==c]
    lodstop <- max(snps$P[snps$ispeak & snps$choose==c]) - loddrop
    snps$minps[snps$choose==c] <- max(c(peakps-maxdist, snps$ps[snps$choose==c & snps$P < lodstop & snps$ps < peakps]))
    snps$maxps[snps$choose==c] <- min(c(peakps+maxdist, snps$ps[snps$choose==c & snps$P < lodstop & snps$ps > peakps]))
  }
  return(snps)
  }else{
  csum <-
    snps %>% group_by(choose) %>% dplyr::summarize(maxps = max(ps), minps = min(ps)) %>% ungroup()
  snps %>% left_join(csum, by = "choose") %>% mutate(maxps = pmin(maxps, ps +
                                                                    maxdist),
                                                     minps = pmax(minps, ps - maxdist))
  }
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
  expp <- ext_peak_sing(pp$gwas, maxdist = args$ldpeakdist)
  write_inrich_phenotype(expp[expp$ispeak == T, ], args$plotdir, phenos[i])
  if (! args$external_inrich){
  run_inrich(
    args$plotdir,
    phenos[i],
    exec = args$inrich,
    i = args$inrich_i,
    j = args$inrich_j
  )
  }
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
    t2 <- ungroup(summarize(group_by(filter(affgen, gene_biotype=="protein_coding"), rs), protein_coding_genes = n()))
    ngene_tbl <- left_join(expp[expp$ispeak == T, ], t2, by = "rs")
    write_csv(ngene_tbl,
              path = paste0(args$plotdir,
                            "/intervals_for_phenotype_",
                            i,
                            ".csv"))
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
  expp <- ext_peak_sing(pp, maxdist = args$ldpeakdist)
  write_inrich_phenotype(expp[expp$ispeak == T, ], args$plotdir, n)
  if (! args$external_inrich){
  run_inrich(
    args$plotdir,
    n,
    exec = args$inrich,
    i = args$inrich_i,
    j = args$inrich_j
  )
  }
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
    # Get the number of protein-coding genes
    t2 <- ungroup(summarize(group_by(filter(affgen, gene_biotype=="protein_coding"), rs), protein_coding_genes = n()))
    ngene_tbl <- left_join(expp[expp$ispeak == T, ], t2, by = "rs")
    if (n == "All Phenotypes"){
      ngene_tbl$groups <- ""
      # Add a column with ';' separated phenotype groups names that overlap the peak
      for (n1 in names(grpwas)){
        if (n1 != "All Phenotypes" && n1 != "All Variance Phenotypes" && n1 != "All Mean Phenotypes"){
          ot <- grpwas[[n1]]
          ot$chr <- as.character(ot$chr)
          ot$chr[ot$chr == "20"] <- "X"
          j1 <- left_join(ngene_tbl, ot[ot$P >= args$pvalthr, ], by = "chr")
          # Filter to where ps is in the minps-maxps range
          j1 <- filter(j1, ps.y >= minps & ps.y <= maxps)
          if (nrow(j1) > 0) {
            print(n1)
            print(j1)
            ngene_tbl$groups[ngene_tbl$rs %in% j1$rs.x] <-
              sapply(ngene_tbl$groups[ngene_tbl$rs %in% j1$rs.x], function(x)
                if (x == "")
                  n1
                else
                  paste(x, n1, sep = ";"))
          }
        }
      }
    }
    write_csv(ngene_tbl, path = paste0(
      args$plotdir,
      "/intervals_for_phenotype_Group_",
      gsub(" ", "_", n),
      ".csv"
    ))
  }
}

for (k in unique(kk$cluster)) {
  write_csv(
    data.frame(genes = clustgene[[k]]),
    path = paste0(args$plotdir, "/genes_for_cluster_", k, ".csv")
  )
  # Run INRICH
  write_inrich_phenotype(clusterpeaks %>% filter(cluster == k),
                         args$plotdir,
                         paste0("cluster_", k))
  if (! args$external_inrich){
    run_inrich(
    args$plotdir,
    paste0("cluster_", k),
    exec = args$inrich,
    i = args$inrich_i,
    j = args$inrich_j
  )}


}

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

# Plot the effectplots for all the peaks
if (args$ploteffect){
  plot_effect(args$outdir, args$plotdir, allpeaks, pnames, fullw, fheight, ffam)
}
