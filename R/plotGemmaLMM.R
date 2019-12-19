#' Plot the GWAS results as a Manhattan plot and highlight specific genes
#'
#' @param results_file A GEMMA results file
#' @param name The title
#' @param metasoft set TRUE if the input is metaSOFT output
#' @param pyLMM TRUE if the input is pyLMM output with rs ID in the first column SNP_ID
#' @param annotations If metasoft is TRUE then annoattions file should be given
#' @param namethr Print gene name above this threshold
#' @param redthr Red points above this thr
#' @param diff A file with results to be subtracted from the first file. Must be in the same format, only implemented for GEMMA
#' @param genotypes The genotypes of the input strains to compute correlation. If given (as data.frame with row.names) every peak will be colored
#' @param maxdist maximal distance between peak and related SNPs
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
plot_gemma_lmm <- function(results_file, name="GWAS results", metasoft=FALSE, pyLMM=FALSE, annotations=NULL, namethr=5, redthr=4, diff=NULL, genotypes=NULL, maxdist=1000000) {
  if (metasoft){
    gwas_results <- read_delim(results_file, "\t", col_names = FALSE, skip=1, guess_max = Inf)
    gwas_results <- gwas_results %>% select(rs=X1, p_wald=X9)  # RSID and PVALUE_RE2
    anno <- read_delim(annotations, ",", col_names = c("rs", "ps", "chr"), guess_max = Inf)
    gwas_results <- left_join(gwas_results, anno, by="rs")
  }else if (pyLMM){
    gwas_results <- read_delim(results_file, "\t", col_names = TRUE, guess_max = Inf)
    gwas_results <- gwas_results %>% select(rs=SNP_ID, p_wald=P_VALUE)  # RSID and PVALUE_RE2
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
    if (! is.null(diff)){
      difres <- read_delim(diff, "\t", col_names = TRUE, col_type = cols(
        .default = col_double(),
        chr = col_character(),
        rs = col_character(),
        ps = col_double(),
        n_miss = col_double(),
        allele1 = col_character(),
        allele0 = col_character()))
      jres <- gwas_results %>% inner_join(select(difres, c(rs, p_wald)), by="rs", suffix = c("", ".d"))
      gwas_results <- jres %>% mutate(p_wald = p_wald/p_wald.d) %>% select(-p_wald.d)
    }
  }

  # Remove correlated peaks
  rep_peaks <- function(genotypes, gwas_pvs, rs_thr=0.6, pthr=1e-20, mxd=maxdist){
    tmat <- base::t(genotypes)
    srt_pv <- gwas_pvs %>% select(rs, p_wald) %>% arrange(p_wald) %>% mutate(choose = 0, ispeak=FALSE)
    peaknum = 1
    while (any(srt_pv$choose[srt_pv$p_wald<=pthr] == 0)){
      nr <- which(srt_pv$choose == 0 & srt_pv$p_wald <= pthr)[1]
      rs <- srt_pv$rs[nr]
      # Select other SNPs in its vicinity
      rsnum <- gwas_pvs$rs==rs
      subt <- tmat[,gwas_pvs$rs[gwas_pvs$chr==gwas_pvs$chr[rsnum] & gwas_pvs$ps >= gwas_pvs$ps[rsnum]-mxd &
                                        gwas_pvs$ps <= gwas_pvs$ps[rsnum]+mxd]]
      if (!is.null(dim(subt))){
        cvec <- cor(tmat[,rs], subt)
        rel_rs <- colnames(cvec)[cvec[1,]^2 >= rs_thr]
        srt_pv[srt_pv$rs %in% rel_rs & srt_pv$choose==0, "choose"] = peaknum
      }
      srt_pv[nr, "choose"] = peaknum
      srt_pv[nr, "ispeak"] = TRUE
      peaknum = peaknum + 1
    }
    return(srt_pv %>% select(rs, choose, ispeak))
  }

  #chr     rs      ps      n_miss  allele1 allele0 af      beta_1  beta_2  beta_3  Vbeta_1_1       Vbeta_1_2       Vbeta_1_3       Vbeta_2_2       Vbeta_2_3       Vbeta_3_3       p_lrt
  #"1"     "rs32166183"    3046097 0       "A"     "C"     0.300   4.737279e-02    1.737096e-02    6.561576e-02    1.160875e-03    9.232757e-04    2.029432e-03    1.757942e-03    2.437142e-03    4.390245e-03    5.048649e-01
  # Add peak color if genotypes are supplied
  genesdist = maxdist/2
  if (!is.null(genotypes)){
    #allgeno <- read.csv(genotypes, header = FALSE, row.names = 1)
    #allgeno <- allgeno[, 3:ncol(allgeno)]
    pnums <- rep_peaks(genotypes, gwas_results, pthr=10^-redthr)
    gwas_results <- gwas_results %>% left_join(pnums, by="rs")
    genesdist = 1000
  }else{
    gwas_results <- gwas_results %>% mutate(choose=0, ispeak=FALSE)
  }
  gwas_results <- gwas_results %>% mutate(P=-log10(p_wald))
  ret_gwas <- gwas_results
  gwas_results[gwas_results$chr=="X","chr"] <- 20# gwas_results %>% dplyr::filter(chr=="X") %>% dplyr::mutate(chr=20)
  gwas_results <- gwas_results %>% mutate(chr=as.numeric(chr)) %>% arrange(chr, ps)

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
    mutate(BPcum=ps+tot) %>%

    # Add highlight and annotation information
    #mutate( is_highlight=ifelse(SNP_ID %in% snpsOfInterest, "yes", "no")) %>%

    # Filter SNP to make the plot lighter
    filter(! is.na(chr))
    # Replace chr 20 to X
 # don[don$chr==20, "chr"] = "X"


  # Prepare X axis
  axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  # Get the RSIDs to put names on
  # Prepare text description for each SNP:
  don$text <- paste("SNP: ", don$rs, "\nPosition: ", don$ps, "\nChromosome: ", don$chr, "\nLOD score:", don$P %>% round(2), "\nWhat else do you wanna know", sep="")
  log10P <- don$P
  ymax <- 1.25 * max(log10P, na.rm = TRUE)
  ymin <- 1.25 * min(log10P, na.rm = TRUE)
  chr_label <- axisdf$chr
  chr_label[chr_label==20] = "X"
  # Make the plot
  p <- ggplot2::ggplot(don, aes(x=BPcum, y=P)) +

    # Show all points
    geom_point(aes(color=as.factor(chr * ((!ispeak)+0) + (((ispeak)+0) * (20+choose) ))) , alpha=1, size=1) +
    scale_color_manual(values = c(rep(c("#CCCCCC", "#969696"),10), rep(c("#cc0029", "#00cc4e", "#0022cc", "#aa00cc"), ceiling(max(gwas_results$choose)/4)) )) +

    # custom X axis:
    scale_x_continuous( label = chr_label, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    ylim(ymin,ymax) +
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
  #if (!is.null(genes)){
    #gwas_results <-  left_join(gwas_results, genes, by="rs")
  names_to <- gwas_results %>% group_by(choose) %>% dplyr::summarize(maxps = max(ps), minps = min(ps)) %>% ungroup()
  names_to <- gwas_results %>% left_join(names_to, by="choose")
  toprs <- get_genes(names_to %>% filter(P>namethr, ispeak==TRUE), dist = genesdist) %>%
    filter(!is.na(mgi_symbol), !stringr::str_detect(mgi_symbol, "Rik$"), !stringr::str_detect(mgi_symbol, "^Gm")) %>%
    group_by(mgi_symbol, chr) %>% summarize(rs=rs[which.max(P)]) %>%
    # Select only one gene
    group_by(rs) %>% summarize(mgi_symbol=mgi_symbol[1])
    # Add gene_name to don
  p <- p + ggrepel::geom_text_repel(data = dplyr::filter(don, rs %in% toprs$rs) %>% left_join(toprs, by="rs"),
                                    aes(BPcum, P, label = mgi_symbol), alpha = 0.7)


  return(list(plot=p, gwas=ret_gwas, pwas=don))
}

#' Title
#'
#' @param snps snps table
#' @param chr chromosome to plot
#' @param from_mb from position in mega bp
#' @param to_mb to position in mega bp
#'
#' @return
#' @export
#' @import qtl2
#' @examples
plot_zoom <- function(snps, chr, from_mb, to_mb, width=1){
  # Download the genes data if not found
  if (!file.exists("mouse_genes_mgi.sqlite")){
    download.file("https://ndownloader.figshare.com/files/17609252", "mouse_genes_mgi.sqlite")
  }
  query_genes <- create_gene_query_func("mouse_genes_mgi.sqlite")
  genes_reg <- query_genes(chr, from_mb, to_mb)
}
