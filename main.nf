#!/usr/bin/env nextflow

/*
 * A pipeline for running mousegwas with shuffling and figure generation
 * https://github.com/TheJacksonLaboratory/mousegwas
 *
 * Parameters:
 *   - input: csv input file with phenotypes
 *   - yaml: yaml file that describes the run
 *   - outdir: output directory for all scripts
 *   - clusters: Number of clusters for clustering QTL (default: 7)
 *   - downsample: Number of individuals to sample from each strain (default: 10)
 *   - genotype: either genotyping files (with -g) or --MDA (default)
 *   - shuffles: Number of shuffles for computing empirical p-values (default: 1000)
 *   - shufyaml: A yaml file describing the shuffle run, should include one phenotype
 *   - pvalue: Threshold for calling a QTL (default 0.05)
 *   - addpostp: Additional parameters for postprocess_mv.R
 *   - addgwas: Additional parameters for run_GWAS.R
 */

params.yaml = ""
params.outdir = "."
params.clusters = 7
params.downsample = 10
params.genotype = "--MDA"
params.shuffles = 1000
params.shufyaml = ""
params.pvalue = 0.05
params.addpostp = ""
params.addgwas = ""
params.input = "NO_FILE"
input = file(params.input)
Channel.fromPath(params.yaml).into{yaml; yaml2}
Channel.fromPath(params.shufyaml).set{shufyml}
process GWAS{
  label 'mousegwas'
  label 'long_run'
  publishDir path:params.outdir, mode:'copy'
  input:
    file infile from input
    file yamfile from yaml

  output:
    file "outdir" into outdir, outdir2

  script:
  def instr = infile.name != "NO_FILE" ? "-i $infile" : ''
  """
  Rscript -e 'source(file=system.file("exec/run_GWAS.R", package="mousegwas"))' $instr -y $yamfile --basedir outdir -d ${params.downsample} ${params.genotype} ${params.addgwas}
  """
}

Channel.of(1..params.shuffles).set{shuffs}
process shuffle{
  label 'mousegwas'
  label 'single_cpu'
  input:
    file yaml from shufyml.collect()
    file infile from input
    val shnum from shuffs
  output:
    file "best_${shnum}.txt" into shout
  script:
  def instr = infile.name != "NO_FILE" ? "-i $infile" : ''
  """
  Rscript -e 'source(file=system.file("exec/run_GWAS.R", package="mousegwas"))' $instr -y $yaml --basedir outdir -d ${params.downsample} --nomv ${params.genotype} --shuffle ${params.addgwas} --seed ${shnum}
  cut -f 13 outdir/output/lmm_pheno_*_all_LOCO.assoc.txt |grep -v p_wald |sort -k1g |head -1 > best_${shnum}.txt
  """
}

process collectShf{
  publishDir path:params.outdir, mode:'copy'
  input:
    file allsh from shout.collect()
  output:
    stdout into threshold
    file "sorted_shuffled_pvalues.txt" into shflist

  script:
  int heads = params.shuffles * params.pvalue
  """
  cat $allsh | sort -k1g |tee sorted_shuffled_pvalues.txt| head -n $heads |tail -1
  """
}

process postp{
  label 'mousegwas'
  label 'high_mem'
  publishDir path:params.outdir, mode:'copy'
  input:
    file yml from yaml2
    file outdir from outdir
    val thr from threshold
  output:
    file "pvalue-threshold.txt" into pvalt
    file "postprocess_nomv/*" into outfiles2
    file "postprocess_nomv/genes_coordinates_for_INRICH.txt" into ggCh
    file "postprocess_nomv/SNPs_map_for_INRICH.txt" into snpCh
    file "postprocess_nomv/intervals*INRICH.txt" into interLs
    file "postprocess_nomv/groups*INRICH.txt" into groupsLs
  script:
  """
  outval=$thr
  echo \$outval > pvalue-threshold.txt
  Rscript -e 'source(file=system.file("exec/postprocess_mv.R", package="mousegwas"))' -p postprocess_nomv -c ${params.clusters} --external_inrich -s 10000 --nomv --pvalthr \$(awk -v o=\$outval 'BEGIN{print -(log(o)/log(10))}') --peakcorr 0.2 -o $outdir -i 20 -j 200 -y $yml ${params.addpostp}
  """
}

//Transforms the lists into channels
interLs.flatMap().set{interCh}
groupsLs.flatMap().set{groupsCh}

//Execute each INRICH run separately
process runINRICH{
  publishDir path:params.outdir, mode:'copy'
  label 'mousegwas'
  label 'single_cpu'
  input:
    file groups from groupsCh
    each file(interv) from interCh
    each file(genes) from ggCh
    each file(snps) from snpCh
  output:
    file "*.out.inrich" into inrout

  script:
  """
  inrich  -c -a $interv -m $snps -g $genes -t $groups -o ${interv.baseName}_${groups.baseName} -i 20 -j 200
  """
}

