#' Script for executing GEMMA and reading its output

#' Download gemma if not available, return the executable
#'
#' @param version gemma version to download
#'
#' @return the gemma path
#' @export
#'
#' @examples
get_gemma <- function(basedir, version = "0.98.1"){
  exec <- Sys.which("gemma")[1]
  if (exec == ""){
    system(paste0("curl -L https://github.com/genetics-statistics/GEMMA/releases/download/",version,"/gemma-",version,"-linux-static.gz |zcat - > ",basedir,"/gemma"))
    system(paste0("chmod a+x ",basedir, "/gemma"))
    exec <- paste0(basedir,"/gemma")
  }
  return(exec)
}


#' Call gemma to calculate kinship matrix
#'
#' @param genotypes A data.table with genotypes of all mice
#' @param annot Annotations of markers, including chr column
#' @param exec gemma executable
#' @param chrname chr for LOCO
#' @param basedir base directory to write files to. Will override existing files
#' @param phenofile a file with phenotypes matching the genotypes table
#' @return The kinship file name
#' @export
#'
#' @importFrom data.table merge fwrite
#' @examples
calc_kinship <- function(genotypes, annot, exec, chrname, basedir, phenofile){
  loco_geno <- genotypes[genotypes$rs %in% annot[annot$chr!=chrname,"rs"],]
  print(head(loco_geno))
  # Write the genotypes without the chr to csv file
  locofname <- paste0(basedir, "/genotypes_LOCO_chr_", chrname, ".csv")
  fwrite(loco_geno, locofname, col.names=FALSE, na="NA")
  # Write a dummy phenotypes file
  # Execute kinship calc in gemma
  system(paste0("cd ", basedir, " && ", exec, " -g ", locofname,
                " -p ", phenofile, " -gk 1 -o kinship_loco_", chrname))
  return(paste0(basedir, "/output/kinship_loco_", chrname, ".cXX.txt"))
}


#' Run lmm 2 (LRT) on the data with LOCO
#'
#' @param genotypes the genotypes data.table
#' @param phenotypes the phenotypes data.table
#' @param annot annotations of SNPs file
#' @param covars covariates data.table
#' @param exec gemma executable
#' @param basedir dir to write files to
#'
#' @return The unified output file
#' @export
#'
#' @importFrom data.table merge fwrite setkey
#' @examples
execute_lmm <- function(genotypes, phenotypes, annot, covars, basedir){
  exec <- get_gemma(basedir)
  # Set keys and merge the genotypes and annotations
  setkey(genotypes, rs, physical = FALSE)
#  setkey(annot, rs, physical = FALSE)
#  loco_geno <- merge(genotypes, annot, by="rs", all.x=TRUE, all.y=FALSE)

  # Write files to disk
  dir.create(basedir, recursive = TRUE)
  phenofile <- paste0(basedir, "/phenotypes.csv")
  fwrite(phenotypes, phenofile, col.names = FALSE, na = "NA", sep=",")
  anotfile <- paste0(basedir, "/annotations.csv")
  fwrite(annot, anotfile, col.names = FALSE, na = "NA", sep=",")
  covarfile <- paste0(basedir, "/covariates.csv")
  fwrite(covars, covarfile, col.names = FALSE, na="NA")
  genofile <- paste0(basedir, "/all_genotypes.csv")
  fwrite(genotypes, genofile, col.names = FALSE, na = "NA")
  # Compute kinship to each chromosome and run gemma with loco
  for (chr in unique(annot$chr)){
    ksfile <- calc_kinship(genotypes, annot, exec, chr, basedir, phenofile)
    geno_sfile <- paste0(basedir, "/genotypes_only_chr_", chr, ".csv")
    fwrite(genotypes[genotypes$rs %in% annot[annot$chr==chrname,"rs"],], geno_sfile, col.names=FALSE, na="NA")
    system(paste0("cd ", basedir, " && ", exec, " -lmm 2 -g ", geno_sfile,
                  " -p ", phenofile, " -a ", anotfile,
                  " -c ", covarfile, " -k ", ksfile, " -o lmm_", chr,
                  " -n ", do.call(paste, c(as.list(1:dim(phenotypes)[2], sep=" ")))))
  }
  # Concatenate all loco files into a single output file
  system(paste0("cd ", basedir,
         " && cat output/lmm*.assoc.txt |head -1 > output/all_lmm_associations.assoc.txt",
         " && cat output/lmm*.assoc.txt |grep -v allele >> output/all_lmm_associations.assoc.txt"))
  return(paste0(basedir, "/output/all_lmm_associations.assoc.txt"))
}
