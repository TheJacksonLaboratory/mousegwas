#!/usr/bin/env Rscript
#
# Collect INRICH results and plot a heatmap with terms on the x axis, phenotypes on the y axis and
# P-values as color
library(dplyr)
library(tidyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(yaml)
library(RColorBrewer)
library(argparse)
library(tibble)
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--files", "-f", nargs='+',
                    help = "a list of INRICH results files")
parser$add_argument("--outdir", "-o", default = ".",
                    help="Plots output dir")
parser$add_argument("-y", "--yaml",
                    help="Yaml file which includes the groups under phenotypes")
parser$add_argument("--clusters",
                    '-c',
                    default = 7,
                    type = "integer",
                    help = "Number of peaks clusters that was used in postprocess")
parser$add_argument("--minpval", "-p",
                    type="double",
                    default=0.5,
                    help="Minimal p-value to display the category")
parser$add_argument("--meanvariance", action="store_true", default=FALSE,
                    help="Plot the PVE of mean and variance phenotypes in different plots")

args <- parser$parse_args()
dir.create(args$outdir, recursive = TRUE)
parse.file <- function(fname){
  con <- pipe(paste0('grep "_O1" ', fname))
  a <- try(read.table(con, sep="\t", header = TRUE), silent = TRUE)
  if (class(a) == "try-error"){
    a = NULL
  }
  return(a)
}
yamin <- yaml.load_file(args$yaml)
grpcol <- RColorBrewer::brewer.pal(8, "Accent")
# Read the phenotypes from the yaml file to get group name
pnames <- data.frame(PaperName = character(0), Group = character(0), stringsAsFactors = FALSE)
pheno_names <- c()
for (n in names(yamin$phenotypes)) {
  pheno_names <- c(pheno_names, yamin$phenotypes[[n]]$papername)
  pnames <-
    rbind(
      pnames,
      data.frame(
        Group = if (length(yamin$phenotypes[[n]]$group)) yamin$phenotypes[[n]]$group else "NoGroup",
        PaperName = yamin$phenotypes[[n]]$papername, stringsAsFactors = FALSE
      )
    )
}
pnames <- rbind(pnames, data.frame(Group = "General", PaperName = "All Phenotypes"))
pheno_names <- c(pheno_names, "All_Phenotypes")
if (args$clusters > 1){
  for (c in 1:args$clusters){
    pnames <- rbind(pnames, data.frame(Group = "General", PaperName = paste0("Cluster ", c)))
    pheno_names <- c(pheno_names, paste0("cluster_", c))
  }
}

groupsOrder = if (length(yamin$groups)) c(yamin$groups, "General") else unique(pnames$Group)
if (args$meanvariance){
  pnames <- rbind(pnames, data.frame(Group = c("General", "General"), PaperName = c("All Mean phenotypes", "All Variance phenotypes")))
  pheno_names <- c(pheno_names, "All_Mean_Phenotypes", "All_Variance_Phenotypes")

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
# intervalsStrideLength.20cm_for_INRICH_groups_MP_terms_for_INRICH.out.inrich
# path <- file.path(args$outdir, "intervals*_for_INRICH*INRICH.out.inrich")
allfiles = args$files #Sys.glob(path)
terms = sub(".*groups_(.*)_for_INRICH.*", "\\1",  allfiles, perl=T)
phenotypes = sub(".*intervals(.*)_for_INRICH_groups.*", "\\1",  allfiles, perl=T)
ptbl <- NULL
for (i in 1:length(allfiles)){
  ptblt <- parse.file(allfiles[i])
  if (!is.null(ptblt)){
    ptblt$group <- terms[i]
    ptblt$phenotype <- phenotypes[i]
    ptbl <- rbind(ptbl, ptblt)
  }
}

# ptbl has the results of all the data + phenotype and groups.
# Filter the terms according to minimal p-value
print(table(ptbl$phenotype))

ptbl <- as_tibble(ptbl)
targets <- unique(as.character(ptbl$TARGET[ptbl$PCORR <= args$minpval]))
ptbl <- filter(ptbl, TARGET %in% targets, phenotype %in% rownames(pnames))
ptbl$logpval <- -log10(ptbl$PCORR)
print(table(ptbl$phenotype))
print(tail(ptbl))

# A function to plot all phenotypes x all
plot.group <- function(tbl){
  lgp_range <- seq(0, ceiling(max(tbl$logpval)))
  lbl_range <- 10^(-lgp_range)
  phn <- unique(as.character(tbl$phenotype[tbl$P<args$minpval]))
  tar <- unique(as.character(tbl$TARGET[tbl$P<args$minpval]))
  plt <- tbl %>% select(TARGET, phenotype, logpval) %>% filter(phenotype %in% phn, TARGET %in% tar) %>% pivot_wider(id_cols = phenotype, names_from = TARGET, values_from = logpval) %>% column_to_rownames(var = "phenotype") %>% as.matrix()
  rnames <- structure(pnames[rownames(plt), "PaperName"], names = rownames(plt))
  cvec = pnames[rownames(plt), "color"]
  cvec <- setNames(cvec, rownames(plt))
  ha <- rowAnnotation(Group = rownames(plt), col=list(Group = cvec), show_legend = FALSE)

  h <- Heatmap(plt, row_labels = rnames[rownames(plt)], column_title = "Terms", row_title = "Phenotypes", name = "p-value", heatmap_legend_param = list(at = lgp_range, labels = lbl_range), right_annotation = ha)

  draw(h, padding = unit(c(1.5, 0.1, 0.1, 0.2), "inches"))
}

# Plot the heatmap for each group
for (g in unique(ptbl$group)){
  cairo_pdf(paste0(args$outdir, "/heatmap_", g, ".pdf"), width = 8.25, height = 11.25)
  plot.group(filter(ptbl, group == g))
  dev.off()
}
