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
#' @return
#' @export
#'
#' @examples
plot_effect <- function(basedir, plotdir, rsnames, fullw = 7.25, height=3.45, ffam = "Arial"){
  # Read the genotypes table
  geno <- read.csv(paste0(basedir, "/all_genotypes.csv"), header = F, row.names = 1)
  geno <- geno[rsnames,]
  pheno <- read.csv(paste0(basedir, "/raw_phenotypes.csv"), row.names = 1)
  for (r in rownames(geno)){
    cairo_pdf(paste0(plotdir, "/effect_plot_", r, ".pdf"), width = fullw, height = height, family = ffam, onefile = TRUE)
    for (p in colnames(pheno)){
      df <- data.frame(base::t(geno[r, -1:-2, drop=F]), pheno[,p,drop=F])
      print(ggplot(df, aes_string(r, p, group=r)) + geom_boxplot() + geom_jitter(alpha = 0.5) + theme_bw() +
        scale_x_continuous(
        breaks = c(0,1,2), labels = c(as.character(geno[r,1]),
                   paste0(geno[r,1], "/", geno[r,2]),
                   as.character(geno[r,2]))))
    }
    dev.off()
  }
}
