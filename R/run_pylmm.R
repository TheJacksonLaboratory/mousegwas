#' Run pylmm
#'
#' @param genotypes All genotypes as a data.table, including three columns of rs, major, minor
#' @param phenotypes All phenotypes (n x p)
#' @param covars Covariates or null
#' @param annot SNPs annotations
#' @param basedir wd
#' @param pymm PyLMM executable
#' @param pylmm_kinship pyLMMKinship executable
#' @param loco Use LOCO or not
#'
#' @return the output file
#'
#' @import data.table
#' @export
run_pylmm <- function(genotypes, phenotypes, annot, covars, basedir, pylmm, pylmm_kinship, loco=TRUE){
  # Write the phenotypes
  if (!is.null(covars)){
    phenotypes <- get_residuals(covars, phenotypes)
  }
  phenofile <- paste0(basedir, "/phenotypes.txt")
  fwrite(base::t(phenotypes), phenofile, col.names = FALSE, na = "NA", sep=" ")

  if (loco) {
    for (chrname in unique(annot$chr)){
      ksfile <- calc_pylmm_kinship(genotypes, annot, pylmm_kinship, chrname, basedir)
      geno_sfile <- paste0(basedir, "/genotypes_only_chr_", chrname, ".txt")
      fwrite(genotypes[genotypes$rs %in% annot[annot$chr==chrname,"rs"], -1:-3], geno_sfile, col.names=FALSE, na="NA", sep=" ")
      annot_file <- paste0(basedir, "/SNPs_only_chr_", chrname, ".txt")
      fwrite(genotypes[genotypes$rs %in% annot[annot$chr==chrname,"rs"], 1, drop=F], annot_file, col.names=TRUE, na="NA", sep=" ")
      # Run pylmm on the chromosome
      run_pylmm_exec(pylmm, geno_sfile, annot_file, phenofile, ksfile, ncol(phenotypes), paste0(basedir, "/output_chr_", chrname))
    }
    # Combine all results
    system(paste0("cd ", basedir, " && head -1 output_chr_1_combined.txt > output_all_chrs_combined.txt && cat output_chr_*_combined.txt | grep -iv beta >> output_all_chrs_combined.txt"))
    return(paste0(basedir,"/otuput_all_chrs_combined.txt"))
  }else{
    # Write the genotypes
    ksfile <- calc_pylmm_kinship(genotypes, annot, pylmm_kinship, NULL, basedir)
    genofile <- paste0(basedir, "/all_genotypes.csv")
    fwrite(genotypes[,-1:-3], genofile, col.names = FALSE, na = "NA", sep=" ")
    annot_file <- paste0(basedir, "/SNPs.txt")
    fwrite(genotypes[, 1], annot_file, col.names=FALSE, na="NA", sep=" ")
    run_pylmm_exec(pylmm, genofile, annot_file, phenofile, ksfile, ncol(phenotypes), paste0(basedir, "/output_all_chrs"))
    return(paste0(basedir, "/output_all_chrs_combined.txt"))
  }
}


#' Calculate kinship using PyLMM
#'
#' @param genotypes
#' @param annot
#' @param pylmm_kinship
#' @param chrname
#' @param basedir
#'
#' @return
#' @export
#'
#' @import data.table
#' @examples
calc_pylmm_kinship <- function(genotypes, annot, pylmm_kinship, chrname, basedir){
  loco_geno <- genotypes[genotypes$rs %in% annot[annot$chr!=chrname,"rs"],]
  # Write the genotypes without the chr to csv file
  locofname <- paste0(basedir, "/genotypes_LOCO_chr_", chrname, ".csv")
  fwrite(loco_geno[,-1:-3], locofname, col.names=FALSE, na="NA", sep=" ")
  # Write a dummy phenotypes file
  # Execute kinship calc in gemma
  print(paste0("cd ", basedir, " && ", pylmm_kinship, " --emmaSNP=", locofname,
               " --emmaNumSNPs=", nrow(loco_geno),  " kinship_loco_", chrname, ".txt"))
  system(paste0("cd ", basedir, " && ", pylmm_kinship, " --emmaSNP=", locofname,
                " --emmaNumSNPs=", nrow(loco_geno),  " kinship_loco_", chrname, ".txt"))
  return(paste0(basedir, "/output/kinship_loco_", chrname, ".txt"))
}


#' Run pyLMM for each phenotype and combine using metaSOFT
#'
#' @param pylmm pylmmGWAS.py executable
#' @param geno_sfile the genotypes file
#' @param annot_file the SNPs names file, to be pasted to the results
#' @param phenofile the phenotypes file
#' @param ksfile kinship matrix file
#' @param nphen number of phenotypes
#' @param output_head Name of output files prefix
#'
#' @return The combined output file
#' @export
#'
#' @import data.table
#' @examples
run_pylmm_exec <- function(pylmm, geno_sfile, annot_file, phenofile, ksfile, nphen, output_head){
  # Run each phenotype (1:nphen)
  for (i in 1:nphen){
    system(paste0(pylmm, " --emmaPHENO=", phenofile, " --emmaSNP=", geno_sfile, " --kfile=", ksfile, " -p ", i-1, " ", output_head, "_", i, ".pyLMM"))
  }

  # Run metaSOFT
  if (nphen > 1){
    outfiles <- sapply(1:nphen, function(n) paste0(output_head, "_", n, ".pyLMM"))
    combine_metaSOFT_pylmm(outfiles, paste0(output_head, "_pasted.txt"), paste0(output_head, "_metasoft.txt"))
    system(paste0("paste ", annot_file, " ", output_head, "_metasoft.txt > ", output_head, "_combined.txt"))
  }else{
    system(paste0("paste ", annot_file, " ", output_head, "_1.pyLMM > ", output_head, "_combined.txt"))
  }
  return(paste0(output_head, "_combined.txt"))
}


#' Title
#'
#' @param infiles
#' @param midfile
#' @param outfile
#' @param version
#'
#' @return
#' @export
#'
#' @import data.table
#' @examples
combine_metaSOFT_pylmm <- function(infiles, midfile, outfile, version="2.0.1"){
  # Download Metasoft snd extracts
  # http://genetics.cs.ucla.edu/meta_jemdoc/repository/2.0.1/Metasoft.zip
  hasmeta <- file.exists(paste0("Metasoft.jar"))
  if (!hasmeta){
    system(paste0("curl -L http://genetics.cs.ucla.edu/meta_jemdoc/repository/",version,"/Metasoft.zip > Metasoft.zip"))
    system(paste0("unzip -uo Metasoft.zip"))
  }

  # Read all the input files and write in the desired format
  # chr     rs      ps      n_miss  allele1 allele0 af      beta_1  Vbeta_1_1   p_lrt
  cmass <- fread(paste0(infiles[1]))
  print(head(cmass))
  #SNP_ID  BETA    BETA_SD F_STAT  P_VALUE
  cmass <- cmass[, c("SNP_ID", "BETA", "BETA_SD")]
  for (n in 2:length(infiles)){
    ctmp <- fread(paste0(infiles[n]))[,c("SNP_ID", "BETA", "BETA_SD")]
    cmass <- merge(cmass, ctmp, by="SNP_ID", all=T, suffixes = c("", paste0(".", n)))
  }
  fwrite(cmass, file=midfile, sep = "\t", col.names = FALSE, row.names = FALSE)

  # Run metasoft
  system(paste0("java -jar Metasoft.jar -mvalue -input ",midfile, " -output ", outfile))
}
