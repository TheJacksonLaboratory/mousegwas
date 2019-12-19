# Get the genes from biomaRt and return the intersecting genes

#' Retrieve genes from biomaRt and return the intersecting genes
#'
#' @param snps A list of SNPs with chr, ps columns
#' @param dist Distance of genes from SNPs
#'
#' @return
#' @export
#'
#' @examples
#' @import biomaRt
#' @import tibble
get_genes <- function(snps, dist=1000000){
  # Get the genes from biomaRt
  #library(biomaRt)
  ensembl <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")
  annot <- biomaRt::getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position",
                 "gene_biotype"), mart=ensembl)
  annot$ext_start <- annot$start_position-dist
  annot$ext_end <- annot$end_position+dist
  # Input might have minps/maxps or just ps
  if (!("minps" %in% names(snps) | "maxps" %in% names(snps))){
    snps <- mutate(snps, maxps=ps, minps=ps)
  }
  annot <- as_tibble(annot) %>% filter(!is.na(chromosome_name), !is.na(ext_start), !is.na(ext_end))
  affgene <- left_join(snps, annot, by = c("chr" = "chromosome_name")) %>%
    filter(maxps >= ext_start, minps <= ext_end)
  return(affgene)
}
