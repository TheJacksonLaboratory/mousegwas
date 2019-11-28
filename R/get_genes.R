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
#' @import fuzzyjoin
#' @import tibble
get_genes <- function(snps, dist=1000000){
  # Get the genes from biomaRt
  library(biomaRt)
  library(fuzzyjoin)
  library(tidyverse)
  ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
  annot <- getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position",
                 "gene_biotype"), mart=ensembl)
  annot$ext_start <- annot$start_position-dist
  annot$ext_end <- annot$end_position+dist
  annot <- as_tibble(annot) %>% filter(!is.na(chromosome_name), !is.na(ext_start), !is.na(ext_end))
  affgene <- left_join(snps, annot, by = c("chr" = "chromosome_name")) %>%
    filter(ps >= ext_start, ps <= ext_end)
  #affgene <- fuzzy_left_join(as_tibble(snps), annot, by=c("chr"="chromosome_name", "ps"="ext_start", "ps"="ext_end"), match_fun=c("==", ">=", "<="))
  # There should probably be a better solution but let's do it with a for loop
  #affgene <- data.frame(rs=character(), gene=character())
  #for (s in 1:nrow(snps)){
  #   fg <- annot[(annot$chromosome_name==as.character(snps[s, "chr", drop=T])) & (annot$start_position>snps[s, "ps", drop=T]-dist) & (annot$end_position<snps[s, "ps", drop=T]+dist), "ensembl_gene_id"]
  #   for (f in fg){
  #     affgene <- rbind(affgene, data.frame(rs=snps[s, "rs"], gene=f))
  #   }
  #}
  return(affgene)
}
