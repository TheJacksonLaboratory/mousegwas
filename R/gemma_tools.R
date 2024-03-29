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
  # Write the genotypes without the chr to csv file
  locofname <- paste0(basedir, "/genotypes_LOCO_chr_", chrname, ".csv")
  fwrite(loco_geno, locofname, col.names=FALSE, na="NA")
  # Write a dummy phenotypes file
  # Execute kinship calc in gemma
  system(paste0("cd ", basedir, " && ", exec, " -g ", locofname,
                " -p ", phenofile, " -gk 1 -o kinship_loco_", chrname))
  return(paste0(basedir, "/output/kinship_loco_", chrname, ".cXX.txt"))
}


#' Get resudials of phenotypes
#' Not used any more
#'
#' @param covars covariates matrix
#' @param phenotypes a table of phenotypes
#'
#' @return a table with residuals of phenotypes
#' @export
#'
#' @examples
get_residuals <- function(covars, phenotypes){
  if (is.null(covars)){
    return (phenotypes)
  }
  resids <- NULL
  for (p in names(phenotypes)){
    print(length(as.data.frame(phenotypes)[,p]))
    print(dim(covars))
    lft <- lm(as.data.frame(phenotypes)[,p] ~ covars[,-1,drop=F], na.action=na.exclude)
    print(lft)
    if (is.null(resids)){
      resids <- data.frame(p=resid(lft))
    }else{
      resids <- cbind(resids, data.frame(resid(lft)))
    }
  }
  colnames(resids) <- names(phenotypes)
  return(resids)
}

#' Combine phenotypes using metasoft
#'
#' @param basedir
#' @param infiles a vector of input files, each for one phenotype
#' @param outfile Output file name
#' @param runit Actually run metasoft. Otherwise only paste the files
#' @param version Metasoft version to download
#' @param xargs additional metasoft args
#'
#' @return Write the unified p_lrt into the outfile
#' @import data.table
#' @export
#'
#' @examples
combine_metaSOFT <- function(basedir, infiles, midfile, outfile, runit = FALSE, version="2.0.1", xargs=""){
  # Download Metasoft snd extracts
  # http://genetics.cs.ucla.edu/meta_jemdoc/repository/2.0.1/Metasoft.zip
  if (runit){
    hasmeta <- file.exists(paste0("Metasoft.jar"))
    if (!hasmeta){
      system(paste0("curl -L http://genetics.cs.ucla.edu/meta_jemdoc/repository/",version,"/Metasoft.zip > Metasoft.zip"))
      system(paste0("unzip -uo Metasoft.zip "))
    }
  }
  infiles2 = c()
  for (f in infiles){
    if (file.exists(paste0(basedir, "/output/", f, ".assoc.txt"))){
      infiles2 <- c(infiles2, f)
    }else{
      print(paste0("Warning: Can't find file ", f))
    }
  }
  infiles <- infiles2
  # Read all the input files and write in the desired format
  # chr     rs      ps      n_miss  allele1 allele0 af      beta_1  Vbeta_1_1   p_lrt
  cmass <- fread(paste0(basedir, "/output/", infiles[1], ".assoc.txt"))
  cmass <- cmass[, c("rs", "beta", "se")]
  for (n in 2:length(infiles)){
    ctmp <- fread(paste0(basedir, "/output/", infiles[n], ".assoc.txt"))[,c("rs", "beta", "se")]
    cmass <- merge(cmass, ctmp, by="rs", all=T, suffixes = c("", paste0(".", n)))
  }
  fwrite(cmass, file=midfile, sep = "\t", col.names = FALSE, row.names = FALSE)
  Sys.sleep(1)
  # Run metasoft
  if (runit){
    system(paste0("java -jar Metasoft.jar  -mvalue_prior_sigma 1 -mvalue -input ",midfile, " -output ", outfile, " -log ", outfile, ".log ", xargs))
  }
}

#' Run lmm 2 (LRT) on the data with LOCO
#'
#' @param genotypes the genotypes data.table
#' @param phenotypes the phenotypes data.table
#' @param annot annotations of SNPs file
#' @param covars covariates data.table
#' @param exec gemma executable
#' @param basedir dir to write files to
#' @param groups A list of vectors, each containing phenotypes indices to run as multivariate. The return file will be the first group unless runmetasoft is TRUE
#' @param loco perform LOCO (default FALSE)
#' @param runmetasoft Run metasoft. Return the output file.
#' @param metasoft_args additional metasoft arguments
#'
#' @return The unified output file
#' @export
#'
#' @import data.table
#' @examples
execute_lmm <- function(genotypes, phenotypes, annot, covars, basedir, groups = NULL, loco=TRUE, runmetasoft=FALSE, metasoft_args=""){
  exec <- get_gemma(basedir)
  dir.create(paste0(basedir, "/output"), recursive = TRUE, showWarnings = FALSE)
  # Write files to disk
  anotfile <- paste0(basedir, "/annotations.csv")
  fwrite(annot, anotfile, col.names = FALSE, na = "NA", sep=",")
  covar_flg <- ""
  if (! is.null(covars)){
    covarfile <- paste0(basedir, "/covariates.txt")
    fwrite(covars, covarfile, col.names = FALSE, na="NA", quote = FALSE, sep = "\t")
    covar_flg <- paste0(" -c ", covarfile)
  }
  genofile <- paste0(basedir, "/all_genotypes.csv")
  fwrite(genotypes, genofile, col.names = FALSE, na = "NA")

  # Write the phenotype files
  phenofile <- paste0(basedir,"/phenotypes.csv")
  fwrite(phenotypes, phenofile, col.names = FALSE, sep=",", na="NA")

  # Compute lmm without LOCO
  if (!loco){
    # Compute kinship matrix
    system(paste0("cd ", basedir, " && ", exec, " -g ", genofile,
                  " -p ", phenofile, " -gk 1 -o kinship_all"))
    ksfile <- paste0(basedir, "/output/kinship_all.cXX.txt")
    for (n in 1:dim(phenotypes)[2]){
      system(paste0("cd ", basedir, " && ", exec, " -lmm 4 -g ", genofile,
                      " -p ", phenofile, " -a ", anotfile,
                      covar_flg,
                      " -k ", ksfile, " -o lmm_all_phenotype_", n,
                      " -n ", n))
      }
    outfiles <- sapply(1:dim(phenotypes)[2], function(n) paste0("lmm_all_phenotype_", n))
    # Run multivariate for each group
    for (i in names(groups)){
      ns <- paste(groups[[i]], collapse=" ")
      system(paste0("cd ", basedir, " && ", exec, " -lmm 3 -g ", genofile,
                    " -p ", phenofile, " -a ", anotfile,
                    covar_flg,
                    " -k ", ksfile, " -o lmm_all_phenotypes_", i,
                    " -n ", ns))

    }
    if (length(outfiles)>1){
      combine_metaSOFT(basedir, outfiles, paste0(basedir, "/output/lmm_all_noloco_allpheno.assoc.pasted.txt"),
                       paste0(basedir, "/output/lmm_all_noloco_allpheno.assoc.txt"), runit = runmetasoft, xargs=metasoft_args)
      if (runmetasoft)
        return(paste0(basedir, "/output/lmm_all_allpheno.assoc.txt"))
      else if (!is.null(groups))
        return (paste0(basedir, "/output/lmm_all_phenotypes_", names(groups)[1], ".assoc.txt"))
    }else
      return(paste0(basedir, "/output/lmm_all_phenotype_1.assoc.txt"))
  }else{ # LOCO
    # Compute kinship to each chromosome and run gemma with loco
    for (chrname in unique(annot$chr)){
      # Avoid empty chromosomes (like MT)
      if (length(intersect(genotypes$rs, annot[annot$chr==chrname,"rs"]))==0) next
      ksfile <- calc_kinship(genotypes, annot, exec, chrname, basedir, phenofile)
      geno_sfile <- paste0(basedir, "/genotypes_only_chr_", chrname, ".csv")
      fwrite(genotypes[genotypes$rs %in% annot[annot$chr==chrname,"rs"],], geno_sfile, col.names=FALSE, na="NA")
      # Run each phenotype first
      outfiles <- sapply(1:dim(phenotypes)[2], function(n) paste0(
        "lmm_", chrname, "_pheno_", n))
      for (n in 1:dim(phenotypes)[2]){
        system(paste0("cd ", basedir, " && ", exec, " -lmin 0.01 -lmax 100 -lmm 4 -g ", geno_sfile,
                      " -p ", phenofile, " -a ", anotfile, covar_flg,
                      " -k ", ksfile, " -o ", outfiles[n],
                      " -n ", n))
      }
      # Run on groups
      for (i in names(groups)){
        ns <- paste(groups[[i]], collapse=" ")
        system(paste0("cd ", basedir, " && ", exec, " -lmm 3 -g ", geno_sfile,
                      " -p ", phenofile, " -a ", anotfile,
                      covar_flg,
                      " -k ", ksfile, " -o lmm_",chrname, "_phenotypes_", i,
                      " -n ", ns))
      }
      # If singles combine the results to one file
      if (length(outfiles)>1){
        combine_metaSOFT(basedir, outfiles, paste0(basedir, "/output/lmm_", chrname, "_allpheno.assoc.pasted.txt"),
                         paste0(basedir, "/output/lmm_", chrname, "_allpheno.assoc.txt"),
                         runit = runmetasoft, xargs = metasoft_args)
      }
    }
      # Concatenate all loco files into a single output file
    if (runmetasoft){
      system(paste0("cd ", basedir,
                    " && cat output/lmm*allpheno.assoc.txt |head -1 > output/all_lmm_LOCO_associations.assoc.txt",
                    " && cat output/lmm*allpheno.assoc.txt |grep -i -v beta >> output/all_lmm_LOCO_associations.assoc.txt"))
    }
    if (length(outfiles)>1){
      system(paste0("cd ", basedir,
                    " && cat output/lmm*allpheno.assoc.pasted.txt > output/all_lmm_LOCO_associations.assoc.pasted.txt"))
    }
    for (i in names(groups)){
      system(paste0("cd ", basedir,
                    " && cat output/lmm*_phenotypes_", i, ".assoc.txt |head -1 > output/lmm_phenotypes_", i, "_all_LOCO.assoc.txt",
                    " && cat output/lmm*_phenotypes_", i, ".assoc.txt | grep -v beta >> output/lmm_phenotypes_", i, "_all_LOCO.assoc.txt"))

    }
    for (n in 1:dim(phenotypes)[2]){
      system(paste0("cd ", basedir,
                    " && cat output/lmm*_pheno_", n, ".assoc.txt |head -1 > output/lmm_pheno_", n, "_all_LOCO.assoc.txt",
                    " && cat output/lmm*_pheno_", n, ".assoc.txt | grep -v beta >> output/lmm_pheno_", n, "_all_LOCO.assoc.txt"))
    }
    if (runmetasoft)
      return(paste0(basedir, "/output/all_lmm_LOCO_associations.assoc.txt"))
    else if (length(groups)){
      return(paste0(basedir, "/output/lmm_phenotypes_", names(groups)[1], "_all_LOCO.assoc.txt"))
    }else
      return(paste0(basedir, "/output/lmm_phenotype_1_all_LOCO.assoc.txt"))
  }
}

#' Return the beta of each strain. Strains are similar if identical in SNPs
#'
#' @param strains_genomes A table of mice genotypes
#' @param phenotypes A table of mice phenotypes
#' @param covars Covariates to include in lm
#' @param downsample Maximal number of representatives. If 0 use average
#' @param sex A vector which assign sex to each individual, to be crossed with strain
#' @param strain A vector with strain name. If given use it instead of matching.
#'
#' @return A list with $genotypes and $phenotypes
#'
#' @export
#'
#' @examples
average_strain <- function(strains_genomes, phenotypes, covars, downsample, sex, strain=NULL){
  # Select random rows to compare, saves time
  set.seed(100)
  cret <- NULL
  if (is.null(strain)){
    sgs <- rbind(as.list(c(1,1,1, sex)), strains_genomes, use.names=FALSE)
    grows <- c(1, sample(nrow(sgs), min(nrow(sgs), 1000)))

    # Find similar genomes
    genidx <- match(sgs[grows,], sgs[grows,])
  }else{
    sgs <- paste(c(1,2,3, strain), c(1,2,3, sex))
    genidx <- match(sgs, sgs)
  }
  if (downsample > 0){ # sample from each strain
    miceidx = 1:3 # 1:3 is rs, minor, major
    for (i in unique(genidx[-1:-3])){
      if (sum(genidx==i) > 1){
        miceidx = c(miceidx, sample(which(genidx==i), min(downsample, sum(genidx==i))))
      }else{
        miceidx = c(miceidx, which(genidx==i))
      }
    }
  }else{ # Pick the first one of each genome
    miceidx = which(!duplicated(genidx))
  }
  gret <- strains_genomes[,miceidx, with=F]

  # Compute the phenotypes
  if (downsample == 0){ # average with lm
    phen2 <- cbind(as.data.frame(phenotypes), as.data.frame(covars[,-1, drop=FALSE]) )
    phen2$strain <- factor(genidx[-1:-3])
    pret <- NULL
    for (pn in colnames(phenotypes)){
      lmout <- lm(as.formula(paste0(pn, " ~ 0 + ", do.call(paste, c(as.list(colnames(covars)[-1]), sep="+")), " + strain ")), phen2)
      print(summary(lmout))
      pret <- cbind(pret, lmout$coefficients[grepl("^strain", names(lmout$coefficients))])
    }
    pret <- setNames(as.data.frame(pret), colnames(phenotypes))
  }else{ # Use the miceidx above
    pret <- phenotypes[miceidx[-1:-3]-3,,drop=F]
    cret <- covars[miceidx[-1:-3]-3,,drop=F]
  }
  return(list(genotypes = gret, phenotypes = pret, indices = miceidx[-1:-3]-3, covars = cret))
}

#' Extract Vg and Ve from GEMMA log file
#'
#' @param logfile GEMMA log file
#'
#' @return Vg and Ve
#' @export
#'
#' @examples
get_sigmas <- function(logfile){
  ## pve estimate in the null model = 0.0253116
  ## se(pve) in the null model = 0.0525698
  ## vg estimate in the null model = 0.0625506
  ## ve estimate in the null model = 1.63363

  con <- file(logfile, "r")
  vg = NULL
  ve = NULL
  pve = NULL
  pvese = NULL
  while(is.null(ve) | is.null(vg)){
    line <- readLines(con, n=1)
    if (grepl("^## pve estimate", line)) pve <- as.numeric(strsplit(line, " = ")[[1]][2])
    if (grepl("^## se\\(pve\\) ", line)) pvese <- as.numeric(strsplit(line, " = ")[[1]][2])
    if (grepl("^## vg estimate", line)) vg <- as.numeric(strsplit(line, " = ")[[1]][2])
    if (grepl("^## ve estimate", line)) ve <- as.numeric(strsplit(line, " = ")[[1]][2])
  }
  return (list(PVE=pve, PVESE=pvese, Vg=vg, Ve=ve))
}

