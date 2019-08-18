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
#'
#' @param covars covariates matrix
#' @param phenotypes a table of phenotypes
#'
#' @return a table with residuals of phenotypes
#' @export
#'
#' @examples
get_residuals <- function(covars, phenotypes){
  resids <- NULL
  for (p in names(phenotypes)){
    lft <- lm(as.data.frame(phenotypes)[,p] ~ covars)
    if (is.null(resids)){
      resids <- data.frame(p=resid(lft))
    }else{
      resids <- cbind(resids, data.frame(p=resid(lft)))
    }
  }
  return(resids)
}

#' Run lmm 2 (LRT) on the data with LOCO
#'
#' @param genotypes the genotypes data.table
#' @param phenotypes the phenotypes data.table
#' @param annot annotations of SNPs file
#' @param covars covariates data.table
#' @param exec gemma executable
#' @param basedir dir to write files to
#' @param eigens Number of eigenvectors to run gemma on
#' @param loco perform LOCO (default FALSE)
#' @param single run each phenotype separately and then run metaSOFT to combine results
#'
#' @return The unified output file
#' @export
#'
#' @importFrom data.table merge fwrite setkey
#' @examples
execute_lmm <- function(genotypes, phenotypes, annot, covars, basedir, eigens, loco=TRUE, single=TRUE){
  exec <- get_gemma(basedir)
  # Set keys and merge the genotypes and annotations
  #setkey(genotypes, rs, physical = FALSE)
#  setkey(annot, rs, physical = FALSE)
#  loco_geno <- merge(genotypes, annot, by="rs", all.x=TRUE, all.y=FALSE)

  # Convert the phenotypes to residuals and do svd to remove correlation
  residuals <- get_residuals(covars, phenotypes)
  resid_comp <- svd(residuals)

  # Write files to disk
  phenofile <- paste0(basedir, "/phenotypes.csv")
  fwrite(resid_comp$u, phenofile, col.names = FALSE, na = "NA", sep=",")
  anotfile <- paste0(basedir, "/annotations.csv")
  fwrite(annot, anotfile, col.names = FALSE, na = "NA", sep=",")
  covarfile <- paste0(basedir, "/covariates.csv")
  fwrite(covars, covarfile, col.names = FALSE, na="NA")
  genofile <- paste0(basedir, "/all_genotypes.csv")
  fwrite(genotypes, genofile, col.names = FALSE, na = "NA")

  if (single){
    # Write the phenotype files
    for (n in range(dim(phenotypes)[2])){
      fwrite(phenotypes[,n,with=FALSE], paste0(basedir,"/phenotype_",n,".csv"), col.names=FALSE, sep=",")
    }
  }

  # Compute lmm without LOCO
  if (!loco){
    system(paste0("cd ", basedir, " && ", exec, " -g ", genofile,
                  " -p ", phenofile, " -gk 1 -o kinship_all"))
    ksfile <- paste0(basedir, "/output/kinship_all.cXX.txt")
    if (single){
      for (n in 1:dim(phenotypes)[2]){
        pfile <- paste0(basedir,"/phenotype_",n,".csv")
        system(paste0("cd ", basedir, " && ", exec, " -lmm 2 -g ", genofile,
                      " -p ", pfile, " -a ", anotfile,
                      " -c ", covarfile,
                      " -k ", ksfile, " -o lmm_all_phenotype_", n,
                      " -n 1"))
      }

    }else{
      system(paste0("cd ", basedir, " && ", exec, " -lmm 2 -g ", genofile,
                    " -p ", phenofile, " -a ", anotfile,
                    " -c ", covarfile,
                    " -k ", ksfile, " -o lmm_all",
                    " -n ", do.call(paste, c(as.list(1:dim(phenotypes)[2], sep=" ")))))
      return(paste0(basedir,"/output/lmm_all.assoc.txt"))
    }
  }else{
    # Compute kinship to each chromosome and run gemma with loco
    for (chrname in unique(annot$chr)){
      ksfile <- calc_kinship(genotypes, annot, exec, chrname, basedir, phenofile)
      geno_sfile <- paste0(basedir, "/genotypes_only_chr_", chrname, ".csv")
      fwrite(genotypes[genotypes$rs %in% annot[annot$chr==chrname,"rs"],], geno_sfile, col.names=FALSE, na="NA")
     # for (n in range(dim(phenotypes)[2])){
      if (!single){
        pfiles <- c(phenofile)
        outfiles <- c(paste0("lmm_", chrname, "_allpheno"))
        nns <- do.call(paste, c(as.list(1:eigens, sep=" ")))
      }else{
        pfiles <- sapply(1:dim(phenotypes)[2], function(n) paste0(basedir,"/phenotype_",n,".csv"))
        outfiles <- sapply(1:dim(phenotypes)[2], function(n) paste0(
          "lmm_", chrname, "_pheno_", n))
        nns <- "1"
      }

      for (n in 1:length(pfiles)){
        print(paste0("Executing: cd ", basedir, " && ", exec, " -lmin 0.01 -lmax 100 -lmm 2 -g ", geno_sfile,
                     " -p ", pfiles[n], " -a ", anotfile,
                     " -k ", ksfile, " -o ", outfiles[n],
                     " -n ", nns))
        system(paste0("cd ", basedir, " && ", exec, " -lmin 0.01 -lmax 100 -lmm 2 -g ", geno_sfile,
                      " -p ", pfiles[n], " -a ", anotfile,
                      #" -c ", covarfile,
                      " -k ", ksfile, " -o ", outfiles[n],
                      " -n ", nns))
      }
      # If singles combine the results to one file
      if (length(pfiles)>1){
        cmass <- fread(paste0(outfiles[1], ".assoc.txt"))
        for (n in 2:length(outfiles)){
          ctmp <- fread(paste0(outfiles[n], ".assoc.txt"))[,.(rs, p_lrt)]

          cmass <- merge(cmass, ctmp, by="rs", all=T, suffixes("", paste0(".", outfiles[n])))
        }
        cmass[paste0("p_lrt.", outfiles[1]) := p_lrt]
        fwrite(cmass, file=paste0("lmm_", chrname, "_allpheno.assoc.txt"))
      }
    }
    # Concatenate all loco files into a single output file
    system(paste0("cd ", basedir,
           " && cat output/lmm*allpheno.assoc.txt |head -1 > output/all_lmm_associations.assoc.txt",
           " && cat output/lmm*allpheno.assoc.txt |grep -v allele >> output/all_lmm_associations.assoc.txt"))
    return(paste0(basedir, "/output/all_lmm_associations.assoc.txt"))
  }
}
