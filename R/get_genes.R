# Get the genes from biomaRt and return the intersecting genes

#' Retrieve genes from biomaRt and return the intersecting genes
#'
#' @param snps A list of SNPs with chr, ps columns. if null returns genes only
#' @param dist Distance of genes from SNPs
#' @param attempts Maximal number of attempts to try biomaRt
#'
#' @return
#' @export
#'
#' @examples
#' @import biomaRt
#' @import tibble
get_genes <- function(snps = NULL,
                      dist = 1000000,
                      attempts = 5, annot=NULL) {
  # Get the genes from biomaRt
  #library(biomaRt)
  if (is.null(annot)) {
    annot <- NULL
    attn <- 1
    while (attn <= attempts && is.null(annot)) {
      attn <- attn + 1
      ensembl <- NULL
      try(ensembl <-
            biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = 'useast.ensembl.org'))
      if (is.null(ensembl))
        try(ensembl <-
              biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
      try(annot <-
            biomaRt::getBM(
              c(
                "ensembl_gene_id",
                "mgi_symbol",
                "chromosome_name",
                "strand",
                "start_position",
                "end_position",
                "gene_biotype",
              ),
              mart = ensembl
            ))
    }
  }
  if (!is.null(snps)) {
    annot$ext_start <- annot$start_position - dist
    annot$ext_end <- annot$end_position + dist
    # Input might have minps/maxps or just ps
    if (!("minps" %in% names(snps) | "maxps" %in% names(snps))) {
      snps <- mutate(snps, maxps = ps, minps = ps)
    }
    annot <-
      as_tibble(annot) %>% filter(!is.na(chromosome_name),
                                  !is.na(ext_start),
                                  !is.na(ext_end))
    affgene <-
      left_join(snps, annot, by = c("chr" = "chromosome_name")) %>%
      filter(maxps >= ext_start, minps <= ext_end)
    return(affgene)
  } else{
    return(annot)
  }
}
