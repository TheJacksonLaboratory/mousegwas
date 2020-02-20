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
  genes$kdesc <- paste0("KEGG map", genes$kpath)
  write_delim(
    genes %>% select(ensembl_gene_id, kpath, kdesc) %>% filter(kpath != ""),
    path = paste0(basedir, "/KEGG_pathway_link_for_INRICH.txt"),
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
                " -o ", name))
}
