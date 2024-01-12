#' Identify the differentially methylated regions using method comb-p.
#' This function identifies the differentially methylated regions by using comb-p method. Users should provide corresponding
#' file and parameters, including the methylation beta matrix file, the phenomenon file with sample information, the type of platform, and the threshold for minimum number of CpGs in a DMR.
#' @param beta The methylation beta values matrix of samples from two groups.
#' @param pd The phenomenon data file with sample information (With Sample group column named Sample_Group).
#' @param minProbes The threshold for minimum number of probes in a DMR.
#' @param arraytype The type of methylation array, including 450K and EPIC array.
#' @return The DMR set detected by comb-p.
#' @export
#' @examples
#' beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
#' pd=read.csv(system.file("extdata","pd.csv",package = 'DMRIntTk'))
#' combp=DMRInt_combp(beta = beta, pd = pd, arraytype = "450K", minProbes = 3)
#'
#'
DMRInt_combp<-
  function(beta, pd, arraytype=c("450K","EPIC"), minProbes){
    #Check the installation of ENmix and other related packages
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
  #Load the arraytype, phenomenon column and other related information
    if(arraytype == "450K"){
      information <- DMRIntTk::data450
    }else{
      information <- DMRIntTk::data850
    }
    pheno<-pd$Sample_Group
    pheno[pheno=="Tumor"]<-0
    pheno[pheno=="Normal"]<-1
    pheno <- as.matrix(pheno)
    results<-cpg.assoc(beta,pheno,large.data=FALSE)
    minProbes=minProbes
    probcg <- rownames(beta)

    #Choose the CpGs obtained in the beta matrix
    cpg <- information[which(information$Name %in% probcg), ]

    combp_beta <- as.data.frame(matrix(NA, nrow(cpg), 5))
    colnames(combp_beta) <- c("chr", "start", "end", "p", "probe")
    combp_beta$chr <- cpg$CHR
    combp_beta$start <- cpg$MAPINFO
    combp_beta$end <- cpg$MAPINFO + 1
    combp_beta$p <- results$results$P.value
    combp_beta$probe <- cpg$Name
    # The following step generates a .*csv file of DMR result file locally,
    #and in order to pass the DMR result back to the user as the result of the function,
    #we need to read the result file.

    combp<-combp(combp_beta,seed = 0.05,region_plot = F,mht_plot = F)
    combp<-read.csv("resu_combp.csv")
    combp<-combp[combp$fdr<0.05&combp$nprobe>minProbes,]
    combp$methodname="combp"
    combp$chr=paste0("chr",combp$chr)
    combp
  }
