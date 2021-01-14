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


#' Title
#'
#' @param genes
#' @param basedir
#'
#' @return
#' @export
#' @import reshape
#'
#' @examples
write_inrich_expression <- function(basedir, genes){
  grpf <- system.file("extdata", "groups_from_41598_2016_BFsrep19274_MOESM4_ESM.txt", package = "mousegwas")
  grpp <- reshape::melt(read.delim(grpf, header=FALSE, sep="\t", stringsAsFactors = FALSE), id.vars = "V1")[,c(1,3)]
  colnames(grpp) <- c("Group", "mgi_symbol")
  grpp <- grpp[grpp$mgi_symbol != "",]
  grpp$grname <- sapply(grpp$Group, function(x) paste(unlist(strsplit(x, " "))[1:2], sep = "_", collapse = "_"))
  # Translate genes mgi_symblo to ensembl
  mgi_ens <- unique(genes[,c("ensembl_gene_id", "mgi_symbol")])
  grpp <- merge(grpp, mgi_ens, by = "mgi_symbol")
  write.table(
    grpp[, c("ensembl_gene_id", "grname", "Group")],
    file = paste0(basedir, "/groups_brain_expression_link_for_INRICH.txt"),
    sep = "\t",
    col.names = F,
    row.names = F,
    quote = 3
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
    unique(genes %>% dplyr::filter(gene_biotype == "protein_coding") %>% dplyr::select(
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
    path = paste0(basedir, "/groups_KEGG_pathway_link_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
  # Read GO mapping and write it too
  gotbl <- genes %>% dplyr::select(ensembl_gene_id, go_id)

  # Get the descriptions
  gotbl <- gotbl[gotbl$go_id!="",]
  dectbl <-
    tibble(
      go_id = unique(gotbl$go_id),
      ont = sapply(unique(gotbl$go_id), function(x) {
        g = GO.db::GOTERM[[x]]
        if (!is.null(g)) g@Ontology else ""
      }),
      desc = sapply(unique(gotbl$go_id), function(x) {
        g = GO.db::GOTERM[[x]]
        if (!is.null(g))
        paste0(GO.db::GOTERM[[x]]@Term, " (", GO.db::GOTERM[[x]]@Ontology, ")")
        else ""
      })
    )
  gotbl <- left_join(gotbl, dectbl, by = "go_id")
  for (ont in c("CC", "BP", "MF")) {
    write.table(
      as.data.frame(gotbl)[gotbl$ont == ont, c("ensembl_gene_id", "go_id", "desc")],
      file = paste0(basedir, "/groups_GO_", ont, "_terms_link_for_INRICH.txt"),
      sep = "\t",
      col.names = F,
      row.names = F,
      quote = 3
    )
    write_inrich_expression(basedir, genes)
  }

  # Download MP mammalian phenotypes annotations from MGI:
  #system(paste0("curl -L http://www.informatics.jax.org/downloads/reports/MGI_PhenotypicAllele.rpt | awk -F\"\t\" '{split($11, sp, \",\"); for (a in sp) print $10\"\t\"sp[a]}' |grep ENS | sort -k2> tmp1"))
  #system(paste0("curl -L http://www.informatics.jax.org/downloads/reports/MP_EMAPA.rpt | cut -f 1,2 | tr \" \" \"-\" | sort -k1 | join -1 2 -2 1 tmp1 - | awk '{print $2\"\t\"$1\"\t\"$3}' >", basedir, "/MP_terms_for_INRICH.txt && rm tmp1"))
  # Download gene-phenotype relationship from Maayan lab DB
  mgitr <- AnnotationDbi::select(org.Mm.eg.db, columns = c("MGI"), keys = genes$ensembl_gene_id, keytype="ENSEMBL")
  mgitr$MGI <- gsub("^MGI:", "", mgitr$MGI)
  mgitr <- unique(mgitr)
  maypf <- system.file("extdata", "gene_attribute_edges.txt", package = "mousegwas")
  mayphen <- read.table(maypf, skip = 1, sep = "\t", quote="\"", header=TRUE)
  mp <- left_join(mayphen, mgitr, by=c("MGI.Accession" = "MGI")) %>% dplyr::select(ENSEMBL, MPID, Phenotype)
  write.table(
    as.data.frame(mp),
    file = paste0(basedir, "/groups_MP_terms_for_INRICH.txt"),
    sep = "\t",
    col.names = F,
    row.names = F,
    quote = 3
  )
  custmp <- system.file("extdata", "MPhenotype_MGenotype.csv", package = "mousegwas")
  mphen <- read_csv(custmp, col_names = c("MGI", "name", "allele", "MPID","MPname"))
  mp2 <- left_join(mphen, mgitr, by=c("MGI")) %>% dplyr::select(ENSEMBL, MPID, MPname)
  write.table(
    as.data.frame(mp2),
    file = paste0(basedir, "/groups_MPMOTOR_terms_for_INRICH.txt"),
    sep = "\t",
    col.names = F,
    row.names = F,
    quote = 3
  )
  return(unique(genes %>% dplyr::select(-go_id)))
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
        " -c -a intervals",
        name,
        "_for_INRICH.txt",
        " -m SNPs_map_for_INRICH.txt ",
        " -g genes_coordinates_for_INRICH.txt",
        " -t groups_KEGG_pathway_link_for_INRICH.txt",
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
        " -c -a intervals",
        name,
        "_for_INRICH.txt",
        " -m SNPs_map_for_INRICH.txt ",
        " -g genes_coordinates_for_INRICH.txt",
        " -t groups_MP_terms_for_INRICH.txt",
        " -o ",
        name,
        "_MP_terms",
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
        " -c -a intervals",
        name,
        "_for_INRICH.txt",
        " -m SNPs_map_for_INRICH.txt ",
        " -g genes_coordinates_for_INRICH.txt",
        " -t groups_MPMOTOR_terms_for_INRICH.txt",
        " -o ",
        name,
        "_MPMOTOR_terms",
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
        " -c -a intervals",
        name,
        "_for_INRICH.txt",
        " -m SNPs_map_for_INRICH.txt ",
        " -g genes_coordinates_for_INRICH.txt",
        " -t groups_brain_expression_link_for_INRICH.txt",
        " -o ",
        name,
        "_brain_expression",
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
          " -c -a intervals",
          name,
          "_for_INRICH.txt",
          " -m SNPs_map_for_INRICH.txt ",
          " -g genes_coordinates_for_INRICH.txt",
          " -t groups_GO_",
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
