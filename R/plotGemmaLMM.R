#' Plot the GWAS results as a Manhattan plot and highlight specific genes
#'
#' @param results_file A GEMMA results file
#' @param genes A tibble with the columns rs and gene_name linking genes to SNPs
#' @param name The title
#' @param metasoft set TRUE if the input is metaSOFT output
#' @param annotations If metasoft is TRUE then annoattions file should be given
#'
#' @return A plot
#' @export
#'
#' @import dplyr
#' @import readr
#' @import ggplot2
#' @importFrom magrittr `%>%`
#' @importFrom readr read_delim
#' @examples
plot_gemma_lmm <- function(results_file, genes=NULL, name="GWAS results", metasoft=FALSE, annotations=NULL) {
  if (metasoft){
    gwas_results <- read_delim(results_file, "\t", col_names = FALSE, skip=1)
    gwas_results <- gwas_results %>% select(rs=X1, p_wald=X9)  # RSID and PVALUE_RE2
    anno <- read_delim(annotations, ",", col_names = c("rs", "ps", "chr"), guess_max = Inf)
    gwas_results <- left_join(gwas_results, anno, by="rs")
  }else{
    gwas_results <- read_delim(results_file, "\t", col_names = TRUE, col_type = cols(
      .default = col_double(),
      chr = col_character(),
      rs = col_character(),
      ps = col_double(),
      n_miss = col_double(),
      allele1 = col_character(),
      allele0 = col_character()))
  }
  #chr     rs      ps      n_miss  allele1 allele0 af      beta_1  beta_2  beta_3  Vbeta_1_1       Vbeta_1_2       Vbeta_1_3       Vbeta_2_2       Vbeta_2_3       Vbeta_3_3       p_lrt
  #"1"     "rs32166183"    3046097 0       "A"     "C"     0.300   4.737279e-02    1.737096e-02    6.561576e-02    1.160875e-03    9.232757e-04    2.029432e-03    1.757942e-03    2.437142e-03    4.390245e-03    5.048649e-01

  gwas_results[gwas_results$chr=="X","chr"] <- 20# gwas_results %>% dplyr::filter(chr=="X") %>% dplyr::mutate(chr=20)
  gwas_results <- gwas_results %>% mutate(chr=as.numeric(chr), P=-log10(p_wald)) %>% arrange(chr, ps)

  #gwasResults <- left_join(gwasResults, snps, by="RSID") %>% left_join(genes, by="RSID")

  #gwasResults$CHROMOSOME[which(gwasResults$CHROMOSOME == "X")] <- "20"
  #gwasResults$CHROMOSOME  <- as.numeric(gwasResults$CHROMOSOME)
  #gwasResults <- gwasResults[order(gwasResults$CHROMOSOME, gwasResults$BASE_POSITION), ]
  #gwasResults$P <- -1 * log10(gwasResults$PVALUE_RE2)

  #nan.ids <- which(is.nan(gwasResults$P))

  #if(length(nan.ids) > 0) {
  #  gwasResults <- gwasResults[-nan.ids,  ]
  #} else {
  #  gwasResults <- gwasResults
  #}

  don <- gwas_results %>%

    # Compute chromosome size
    group_by(chr) %>%
    summarise(chr_len=max(ps)+10000000) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(gwas_results, ., by=c("chr"="chr")) %>%

    # Add a cumulative position of each SNP
    arrange(chr, ps) %>%
    mutate( BPcum=ps+tot) %>%

    # Add highlight and annotation information
    #mutate( is_highlight=ifelse(SNP_ID %in% snpsOfInterest, "yes", "no")) %>%

    # Filter SNP to make the plot lighter
    filter(P>0.5)



  # Prepare X axis
  axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  # Get the RSIDs to put names on
  # Prepare text description for each SNP:
  don$text <- paste("SNP: ", don$rs, "\nPosition: ", don$ps, "\nChromosome: ", don$chr, "\nLOD score:", don$P %>% round(2), "\nWhat else do you wanna know", sep="")
  log10P <- don$P
  ymax <- 1.25 * max(log10P, na.rm = TRUE)
  # Make the plot
  p <- ggplot2::ggplot(don, aes(x=BPcum, y=P)) +

    # Show all points
    geom_point( aes(color=as.factor(chr + 21 * ((P>8)+0))), alpha=0.8, size=0.5) +
    scale_color_manual(values = c(rep(c("grey", "skyblue"),10), rep("red", 20) )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    ylim(0,ymax) +
    xlab(name) +
    ylab("-log(P-value)") +

    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +

    # Custom the theme:
    theme_bw() +
    theme(
      legend.position="none",
      text = element_text(size=20),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  if (!is.null(genes)){
    gwas_results <-  left_join(gwas_results, genes, by="rs")
    toprs <- gwas_results %>% filter(P>8, !is.na(gene_name), !stringr::str_detect(gene_name, "Rik$"), !stringr::str_detect(gene_name, "^Gm")) %>% group_by(gene_name, chr) %>% summarize(rs=rs[which.max(P)]) %>%
      # Select only one gene
      group_by(rs) %>% summarize(gene_name=gene_name[1])
    p <- p + ggrepel::geom_text_repel(data = dplyr::filter(don, rs %in% toprs$rs, gene_name %in% toprs$gene_name),
                                      aes(BPcum, P, label = gene_name), alpha = 0.7)
  }

  print(p)
  return(p)
}
