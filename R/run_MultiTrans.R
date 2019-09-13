#' Download MultiTrans
#'
#' @param weblink Where to download the zip file from
#'
#' @return TRUE/FALSE if download was successful
#' @export
#'
#' @examples
get_multi <- function(multilink="http://genetics.cs.ucla.edu/multiTrans/MultiTrans.zip",
                      slidelink="http://slide.cs.ucla.edu/data/slide.1.0.4.intel64.tar.gz"){
  if (! file.exists("MultiTrans")){
    system(paste0("curl -L ", multilink, " > MultiTrans.zip"))
    system("unzip MultiTrans.zip")
  }
  if (! file.exists("slide.1.0.4.intel64.tar.gz")){
    system(paste0("curl -L ", slidelink, " > slide.1.0.4.intel64.tar.gz"))
    system("tar xvzf slide.1.0.4.intel64.tar.gz")
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
  while(is.null(ve) & is.null(vg)){
    line <- readLines(con, n=1)
    if (grepl("^## vg estimate", line)) vg <- as.numeric(strsplit(line, " = ")[[1]][2])
    if (grepl("^## ve estimate", line)) ve <- as.numeric(strsplit(line, " = ")[[1]][2])
  }
  return (vg, ve)
}

#' Wrap generateR.R
#'
#' @param genotypes The genotypes matrix (just the values)
#' @param kinship The kinship matrix
#' @param vg_med median Vg value
#' @param ve_med median Ve value
#' @param basedir dir to write r.txt to
#'
#' @return
#' @export
#'
#' @examples
run_R <- function(genotypes, kinship, vg_med, ve_med, basedir){
  # Load the source code to use the functions
  source("MultiTrans/generateR.R")

  # This code was modified from generateR.R
  snpNum <- dim(genotypes)[2]
  indiNum <- dim(genotypes)[1]
  I<-matrix(0,nrow=indiNum, ncol=indiNum);
  I[row(I)==col(I)]<-1;
  sigmaM = vg_med*kinship + ve_med*I;
  UX = rotate(genotypes,sigmaM)
  Ur=cor(UX,UX)
  write.table(Ur,paste0(basedir,"/r.txt"),row.names=FALSE, col.names=FALSE, quote=FALSE)
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
  system(exec1)
  exec2 <- paste0("slide.1.0/slide_2run ", basedir, "/slide_prep.out ", basedir, "/slide_run.out ", samps, " ", seed)
  system(exec2)
  exec3 <- paste0("slide.1.0/slide_3sort ", basedir,"/slide_sort.out ", basedir, "/slide_run.out")
  system(exec3)
  # Write a file with `nsnp` pvals
  write.table(matrix(pval, snps, 1),paste0(basedir, "/pval_base.txt"),row.names=FALSE, col.names=FALSE, quote=FALSE)
  exec4 <- paste0("slide.1.0/slide_4correct -t ", basedir, "/slide_sort.out ", basedir, "/pval_base.txt ", basedir, "/pval_thresholds.txt")
  system(exec4)
  return(paste0(basedir, "/pval_thresholds.txt"))
}
