#' Write the intervals for INRICH execution
#'
#' @param snps a data.frame or tibble with the columns chr, minps an maxps describinf intervals
#' @param basedir Where to write the files to
#' @param name Name of the phenotype
#'
#' @return
#' @export
#'
#' @examples
write_inrich_phenotype <- function(snps, basedir, name) {
  # Write the intervals using the minps and maxps columns
  name <- gsub(" ", "_", name)
  write_delim(
    snps[, c("chr", "minps", "maxps")],
    path = paste0(basedir, "/intervals", name, "_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
}

#' Write the SNPs map for INRICH execution
#'
#' @param snps a data.frame or tibble with the columns chr, minps an maxps describinf intervals
#' @param basedir Where to write the files to
#'
#' @return
#' @export
#'
#' @examples
write_inrich_snps <- function(snps, basedir) {
  # Write the snps chr and position
  write_delim(
    snps %>% dplyr::select(chr, bp38),
    path = paste0(basedir, "/SNPs_map_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
}

#' Write the genes, KEGG and GO mapping for INRICH execution
#'
#' @param basedir Where to write the files to
#'
#' @return
#' @export
#'
#' @examples
write_genes_map <- function(basedir) {
  # Write the genes map chr, staret, stop, ID, desc (mgi_symbol)
  genes <- get_genes()
  write_delim(
    unique(genes %>% dplyr::select(
      chromosome_name,
      start_position,
      end_position,
      ensembl_gene_id,
      mgi_symbol
    )),
    path = paste0(basedir, "/genes_coordinates_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
  library(org.Mm.eg.db)
  paths <- AnnotationDbi::select(org.Mm.eg.db, columns = c("PATH"), keys = genes$ensembl_gene_id, keytype="ENSEMBL")
  kgenes <- left_join(genes, paths, by = c("ensembl_gene_id" = "ENSEMBL"))
  #write the KEGG mapping
  kgenes$kdesc <- paste0("KEGG map", kgenes$PATH)
  write_delim(
    kgenes %>% dplyr::select(ensembl_gene_id, PATH, kdesc) %>% filter(PATH != ""),
    path = paste0(basedir, "/KEGG_pathway_link_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
  # Read GO mapping and write it too
  gotbl <- genes %>% dplyr::select(ensembl_gene_id, goslim_goa_accession)

  # Get the descriptions
  gotbl <- gotbl[!is.na(gotbl$go),]
  dectbl <-
    tibble(
      goslim_goa_accession = unique(gotbl$goslim_goa_accession),
      ont = sapply(unique(gotbl$goslim_goa_accession), function(x) {
        GO.db::GOTERM[[x]]@Ontology
      }),
      desc = sapply(unique(gotbl$goslim_goa_accession), function(x) {
        paste0(GO.db::GOTERM[[x]]@Term, " (", GO.db::GOTERM[[x]]@Ontology, ")")
      })
    )
  gotbl <- join_left(gotbl, dectbl, by = "goslim_goa_accession")
  for (ont in c("CC", "BP", "MF")) {
    write.table(
      as.data.frame(gotbl)[gotbl$ont == ont, c("ensembl_gene_id", "goslim_goa_accession", "desc")],
      file = paste0(basedir, "/GO_", ont, "_terms_link_for_INRICH.txt"),
      sep = "\t",
      col.names = F,
      row.names = F,
      quote = 3
    )
  }

  # Download MP mammalian phenotypes annotations from MGI:
  system(paste0("curl -L http://www.informatics.jax.org/downloads/reports/MGI_PhenotypicAllele.rpt | awk -F\"\t\" '{split($11, sp, \",\"); for (a in sp) print $10\"\t\"sp[a]}' |grep ENS | sort -k2> tmp1"))
  system(paste0("curl -L http://www.informatics.jax.org/downloads/reports/MP_EMAPA.rpt | cut -f 1,2 | tr \" \" \"-\" | sort -k1 | join -1 2 -2 1 tmp1 - | awk '{print $2\"\t\"$1\"\t\"$3}' >", basedir, "/MP_terms_for_INRICH.txt && rm tmp1"))

  return(unique(genes %>% dplyr::select(-goslim_goa_accession)))
}


#' Rnu KEGG and GO INRICH for written intervals
#'
#' @param basedir Where the files are
#' @param name The phenotype name
#' @param exec INRICH executable file
#' @param i minimal group size
#' @param j maximal group size
#'
#' @return
#' @export
#'
#' @examples
run_inrich <-
  function(basedir,
           name,
           exec = "inrich",
           i = 5,
           j = 200) {
    # Remove spaces from name
    name <- gsub(" ", "_", name)
    system(
      paste0(
        "cd ",
        basedir,
        " && ",
        exec,
        " -a intervals",
        name,
        "_for_INRICH.txt",
        " -m SNPs_map_for_INRICH.txt ",
        " -g genes_coordinates_for_INRICH.txt",
        " -t KEGG_pathway_link_for_INRICH.txt",
        " -o ",
        name,
        "_KEGG_pathways",
        " -i ",
        i,
        " -j ",
        j
      )
    )
    system(
      paste0(
        "cd ",
        basedir,
        " && ",
        exec,
        " -a intervals",
        name,
        "_for_INRICH.txt",
        " -m SNPs_map_for_INRICH.txt ",
        " -g genes_coordinates_for_INRICH.txt",
        " -t MP_terms_for_INRICH.txt",
        " -o ",
        name,
        "_MP_terms",
        " -i ",
        i,
        " -j ",
        j
      )
    )
    for (ont in c("CC", "BP", "MF")) {
      system(
        paste0(
          "cd ",
          basedir,
          " && ",
          exec,
          " -a intervals",
          name,
          "_for_INRICH.txt",
          " -m SNPs_map_for_INRICH.txt ",
          " -g genes_coordinates_for_INRICH.txt",
          " -t GO_",
          ont,
          "_terms_link_for_INRICH.txt",
          " -o ",
          name,
          "_GO_",
          ont,
          "_terms",
          " -i ",
          i,
          " -j ",
          j
        )
      )
    }

  }
