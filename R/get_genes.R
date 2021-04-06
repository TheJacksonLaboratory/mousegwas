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
      # Hard-coded version GRCm38 because MDA is with this version.
      try(ensembl <-
            biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="nov2020.archive.ensembl.org"))

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
                "go_id"
              ),
              mart = ensembl,
              useCache = FALSE
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


# Get the genes from biomaRt and return the intersecting genes

#' Retrieve genes from biomaRt and return the intersecting genes
#'
#' @param human_genes A list of human ensembl_gene_id
#'
#' @return a 2-column table with human emsembl_gene_id in the first and mouse in the second
#' @export
#'
#' @examples
#' @import biomaRt
convert_to_mouse <- function(human_genes){
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = biomaRt::getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = human_genes , mart = human, attributesL = c("ensembl_gene_id"), martL = mouse, uniqueRows=T)
  return(genesV2)
}
