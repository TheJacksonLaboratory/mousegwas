#' Plot effectplots for each important SNP
#'
#' @param basedir Where the run_GWAS output lies
#' @param plotdir Where the output pdfs should be written to
#' @param rsnames The list of SNPs to be plotted
#' @param fullw pdf width
#' @param height pdf height
#' @param ffam font family
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import cowplot
#' @return
#' @export
#'
#' @examples
plot_effect <- function(basedir, plotdir, rsnames, pnames, fullw = 7.25, height=11-1.25, ffam = "Arial"){
  # Read the genotypes table
  geno <- read.csv(paste0(basedir, "/all_genotypes.csv"), header = F, row.names = 1, check.names = FALSE)
  geno <- geno[rsnames,]
  pheno <- read.csv(paste0(basedir, "/raw_phenotypes.csv"), row.names = 1)
  colnames(pheno) <- pnames$PaperName[match(colnames(pheno),rownames(pnames))]
  grows <- c(sample(nrow(geno), min(nrow(geno), 1000)))
  strains <- match(geno[grows, -1:-2], geno[grows, -1:-2])
  for (r in rownames(geno)){
    cairo_pdf(paste0(plotdir, "/effect_plot_", r, ".pdf"), width = fullw, height = height, family = ffam, onefile = TRUE)
    pvec <- list()
    i=1
    for (p in colnames(pheno)){
      print(p)
      print(pnames)
      print(head(pheno))
      df <- data.frame(base::t(geno[r, -1:-2, drop=F]), pheno[,p,drop=F], strains = factor(strains))
      df[,r] <- as.factor(df[,r])
      pvec[[i]] <- ggplot(df, aes_string(x=r, y=p, group=r)) + geom_violin(aes_string(fill = r), scale="area") +
              geom_boxplot(width=0.1) +
              scale_fill_brewer(palette = "Pastel1") +
              #geom_jitter(alpha = 0.5, height=0, width=0.1) + theme_bw() +
        scale_x_discrete(
        breaks = c(0,1,2), labels = c(as.character(geno[r,1]),
                   paste0(geno[r,1], "/", geno[r,2]),
                   as.character(geno[r,2]))) + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank())
      i <- i + 1
    }
    print(plot_grid(plotlist = pvec, ncol = 4))
    dev.off()
  }
}
