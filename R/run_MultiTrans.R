#' Download MultiTrans
#'
#' @param weblink Where to download the zip file from
#'
#' @return TRUE/FALSE if download was successful
#' @export
#'
#' @examples
get_multi <- function(multilink="http://genetics.cs.ucla.edu/multiTrans/MultiTrans.zip",
                      slidelinux="http://slide.cs.ucla.edu/data/slide.1.0.4.intel64.tar.gz",
                      slidemac="http://slide.cs.ucla.edu/data/slide.1.0.4.MAC.64bit.tar.gz"){
  if (! file.exists("MultiTrans")){
    system(paste0("curl -L ", multilink, " > MultiTrans.zip"))
    system("unzip MultiTrans.zip")
  }
  if (! file.exists("slide.1.0")){
    # Determine platform
    if (Sys.info()["sysname"] == "Linux"){
      system(paste0("curl -L ", slidelinux, " > slide.1.0.4.intel64.tar.gz"))
      system("tar xvzf slide.1.0.4.intel64.tar.gz")
    }else if (Sys.info()["sysname"] == "Darwin"){
      system(paste0("curl -L ", slidemac, " > slide.1.0.4.MAC.64bit.tar.gz"))
      system("tar xvzf slide.1.0.4.MAC.64bit.tar.gz")
      system("mv slide.1.0.64bit slide.1.0")
    }else{
      print("Can't determine system or not Linux/OSX")
      return(FALSE)
    }
  }
  if (! file.exists("MultiTrans/generateR.R")){
    print("Error: Can't find generateR.R, MultiTrans couldn't be downloaded")
    return (FALSE)
  }
  if (! file.exists("slide.1.0/slide")){
    print("Error: Can't find SLIDE, SLIDE couldn't be downloaded")
    return(FALSE)
  }
  return (TRUE)
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
  ## vg estimate in the null model = 0.00981952
  ## ve estimate in the null model = 0.0178361
  con <- file(logfile, "r")
  vg = NULL
  ve = NULL
  while(is.null(ve) | is.null(vg)){
    line <- readLines(con, n=1)
    if (grepl("^## vg estimate", line)) vg <- as.numeric(strsplit(line, " = ")[[1]][2])
    if (grepl("^## ve estimate", line)) ve <- as.numeric(strsplit(line, " = ")[[1]][2])
  }
  return (list(vg=vg, ve=ve))
}

#' Wrap generateR.R
#'
#' @param genotypes_file The genotypes matrix (just the values)
#' @param kinship_file The kinship matrix
#' @param vg_med median Vg value
#' @param ve_med median Ve value
#' @param basedir dir to write r.txt to
#'
#' @return
#' @export
#'
#' @examples
run_R <- function(genotypes, kinship_file, vg_med, ve_med, basedir){
  # Write the genotypes, transposed
  write.table(base::t(genotypes), file=paste0(basedir, "/gmatrix.txt"), col.names = F, row.names = F, quote = F)
  # Load the source code to use the functions
  write.table(data.frame(vg=c(vg_med), ve=c(ve_med)), file = paste0(basedir, "/VC.txt") ,col.names = FALSE, row.names = FALSE, quote=FALSE)
  #"Usage: R CMD BATCH --args -Xpath= -Kpath= -VGpath= -outputPath= generateR.R generateR.log
  system(paste0("R CMD BATCH --args -Xpath=",basedir, "/gmatrix.txt -Kpath=", kinship_file, " -VCpath=", basedir, "/VC.txt -outputPath=", basedir, " MultiTrans/generateR.R" ))
}

#' Run java generateC
#'
#' @param basedir where to reasd and write
#' @param windowSize default to 1000
#'
#' @return
#' @export
#'
#' @examples
run_C <- function(basedir, windowSize=1000){
  exec <- paste0("java -jar MultiTrans/generateC/generateC.jar ", windowSize, " ", basedir, "/r.txt ", basedir, "/c.txt")
  system(exec)
}

#' Run SLIDE
#'
#' @param basedir Where to read and write the files
#' @param nsnp number of SNPs
#' @param pval desired pvalue (default 0.05)
#' @param windowSize same as in run_C
#' @param samps number of random sampling
#' @param seed rand seed
#'
#' @return output thresholds file
#' @export
#'
#' @examples
run_slide <- function(basedir, nsnp, pval=0.05, windowSize=1000, samps=10000000, seed=100){
  exec1 <- paste0("slide.1.0/slide_1prep -C ", basedir, "/c.txt ", windowSize, " ", basedir, "/slide_prep.out")
  print(exec1)
  system(exec1)
  exec2 <- paste0("slide.1.0/slide_2run ", basedir, "/slide_prep.out ", basedir, "/slide_run.out ", sprintf("%d", samps), " ", seed)
  print(exec2)
  system(exec2)
  exec3 <- paste0("slide.1.0/slide_3sort ", basedir,"/slide_sort.out ", basedir, "/slide_run.out")
  print(exec3)
  system(exec3)
  # Write a file with `nsnp` pvals
  write.table(matrix(pval, nsnp, 1),paste0(basedir, "/pval_base.txt"),row.names=FALSE, col.names=FALSE, quote=FALSE)
  exec4 <- paste0("slide.1.0/slide_4correct -t ", basedir, "/slide_sort.out ", basedir, "/pval_base.txt ", basedir, "/pval_thresholds.txt")
  print(exec4)
  system(exec4)
  return(paste0(basedir, "/pval_thresholds.txt"))
}
