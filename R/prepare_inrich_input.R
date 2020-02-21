write_inrich_phenotype <- function(snps, basedir, name) {
  # Write the intervals using the minps and maxps columns
  write_delim(
    snps %>% select(chr, minps, maxps),
    path = paste0(basedir, "/intervals", name, "_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
}

write_inrich_snps <- function(snps, basedir) {
  # Write the snps chr and position
  write_delim(
    snps %>% select(chr, bp38),
    path = paste0(basedir, "/SNPs_map_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
}

write_genes_map <- function(basedir) {
  # Write the genes map chr, staret, stop, ID, desc (mgi_symbol)
  genes <- get_genes()
  write_delim(
    genes %>% select(
      chromosome_name,
      start_position,
      end_position,
      ensembl_gene_id,
      mgi_symbol
    ),
    path = paste0(basedir, "/genes_coordinates_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
  g <- genes
  genes$kpath <- gsub("\\+.*", "", genes$kegg_enzyme)
  #write the KEGG mapping
  genes$kdesc <- paste0("KEGG map", genes$kpath)
  write_delim(
    genes %>% dplyr::select(ensembl_gene_id, kpath, kdesc) %>% filter(kpath != ""),
    path = paste0(basedir, "/KEGG_pathway_link_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
  # Read GO mapping and write it too
  gotrm <- goseq::getgo(genes$ensembl_gene_id, "mm10", "ensGene")
  gotbl <-
    data.table(gene = character(0),
           go = character(0))
  #gtr <- data.table(go=character(0), desc=character(0))
  #setkey(gtr, go)
  #goterm <- GO.db::GOTERM
  for (i in 1:length(gotrm)) {
    print(i)
    for (g in gotrm[[i]]) {
   #   desc <- gtr[g, desc]
  #    if (length(desc)==0) {
  #      desc = paste0(goterm[[g]]@Term, " (", goterm[[g]]@Ontology, ")")
  #      gtr <-  rbindlist(list(gtr, list(g, desc)))
  #      setkey(gtr, go)
  #    }
      desc=g
    #  gotbl <- add_row(gotbl,
      gotbl <- rbindlist(list(gotbl, list(
                       gene = names(gotrm)[i],
                       go = g)))
    }
  }
  # Get the descriptions
  dectbl <-
    data.table(go = unique(gotbl$go),
               desc = sapply(unique(gotbl$go), function(x) {
                 paste0(GO.db::GOTERM[[x]]@Term, " (", GO.db::GOTERM[[x]]@Ontology, ")")
               }))
  gotbl <- merge.data.table(gotbl, dectbl, by="go")
  write_delim(
    gotbl,
    path = paste0(basedir, "GO_terms_link_for_INRICH.txt"),
    delim = "\t",
    col_names = FALSE
  )
  return(g)
}

run_inrich <- function(basedir, name, exec="inrich"){
  system(paste0("cd ", basedir, " && ", exec,
                " -a intervals", name, "_for_INRICH.txt",
                " -m SNPs_map_for_INRICH.txt ",
                " -g genes_coordinates_for_INRICH.txt",
                " -t KEGG_pathway_link_for_INRICH.txt",
                " -o ", name, "_KEGG_pathways"))
  system(paste0("cd ", basedir, " && ", exec,
                " -a intervals", name, "_for_INRICH.txt",
                " -m SNPs_map_for_INRICH.txt ",
                " -g genes_coordinates_for_INRICH.txt",
                " -t GO_terms_link_for_INRICH.txt",
                " -o ", name, "GO_terms"))
}
