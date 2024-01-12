#' Identify the differentially methylated regions using ipDMR method.
#' This function identifies the differentially methylated regions by using ipDMR method. Users should provide corresponding
#' file and parameters, including the methylation beta matrix file, the phenomenon file with sample information, the type of platform, and the threshold for minimum number of CpGs in a DMR.
#' @param beta The methylation beta values matrix of samples from two groups.
#' @param pd The phenomenon data file with sample information (With Sample group column named Sample_Group).
#' @param arraytype The type of methylation array, including 450K and EPIC array.
#' @param minProbes The threshold for minimum number of probes in a DMR.
#' @return The DMR set detected by ipDMR.
#' @export
#' @examples
#' beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
#' pd=read.csv(system.file("extdata","pd.csv",package = 'DMRIntTk'))
#' ipDMR=DMRInt_ipDMR(beta = beta, pd = pd, arraytype = "450K", minProbes = 3)
#'
DMRInt_ipDMR<-
  function(beta, pd, arraytype=c("450K","EPIC"), minProbes){
  #Check the installation of ENmix and other related packages(combp and ipDMR can
  #be executed by different functions from the same packages)
  if(require("ENmix")){
  print("ENmix package is successfully loaded.")
} else {
  print("ENmix package does not exist, trying to install")
  install.packages("ENmix")
  if(require("ENmix")){
    print("ENmix package is successfully installed, loading now")
  } else {
    stop("Installation failed, please go to Bioconductor to install ENmix package at: https://www.bioconductor.org/packages/release/bioc/html/ENmix.html")
  }
}
if(require("CpGassoc")){
  print("CpGassoc package is successfully loaded.")
} else {
  print("CpGassoc package does not exist, trying to install")
  install.packages("CpGassoc")
  if(require("CpGassoc")){
    print("CpGassoc package is successfully installed, loading now")
  } else {
    stop("Installation failed, please go to Bioconductor to install CpGassoc package at: https://www.bioconductor.org/packages/release/bioc/html/CpGassoc.html")
  }
}
library(ENmix)
library(CpGassoc)
if(arraytype == "450K"){
  information <- DMRIntTk::data450
}else{
  information <- DMRIntTk::data850
}
#The pipeline of ipDMR is the same as the combp.
pheno<-pd$Sample_Group
pheno[pheno=="Tumor"]<-0
pheno[pheno=="Normal"]<-1
pheno <- as.matrix(pheno)
results<-cpg.assoc(beta,pheno,large.data=FALSE)
minProbes=minProbes
probcg <- rownames(beta)
cpg <- information[which(information$Name %in% probcg), ]
ipDMR_beta <- as.data.frame(matrix(NA, nrow(cpg), 5))
colnames(ipDMR_beta) <- c("chr", "start", "end", "p", "probe")
ipDMR_beta$chr <- cpg$CHR
ipDMR_beta$start <- cpg$MAPINFO
ipDMR_beta$end <- cpg$MAPINFO + 1
ipDMR_beta$p <- results$results$P.value
ipDMR_beta$probe <- cpg$Name
ipdmr<-ipdmr(ipDMR_beta,seed = 0.05,region_plot = F ,mht_plot = F)
ipdmr<-read.csv("resu_ipdmr.csv")
ipdmr<-ipdmr[ipdmr$fdr<0.05&ipdmr$nprobe>minProbes,]
ipdmr$methodname="ipDMR"
ipdmr$chr=paste0("chr",ipDMR$chr)
ipdmr
}

