#' Identify the differentially methylated regions using ProbeLasso method.
#' This function identifies the differentially methylated regions by using ProbeLasso method. Users should provide corresponding
#' file and parameters, including the methylation beta matrix file, the phenomenon file with sample information, the type of platform, and the threshold for minimum number of probes in a DMR.
#' @param beta The methylation beta values matrix of samples from two groups.
#' @param pd The phenomenon data with sample information (With Sample group column named Sample_Group).
#' @param arraytype The type of methylation array, including methylation 450K and EPIC array.
#' @param minProbes The threshold for minimum number of probes in a DMR.
#' @return The DMR set detected by ProbeLasso.
#' @export
#' @examples
#' beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
#' pd=read.csv(system.file("extdata","pd.csv",package = 'DMRIntTk'))
#' ProbeLasso=DMRInt_ProbeLasso(beta = beta, pheno=pd$Sample_Group, arraytype = "450K", minProbes = 3)
#'
DMRInt_ProbeLasso<-
  function(beta, pheno, arraytype=c("450K", "EPIC"), minProbes){
    #Check the installation of ChAMP(Probelasso and bumphunter are executed with different functions of the CHAMP package)
    if(require("ChAMP")){
      print("ChAMP package is successfully loaded.")
    } else {
      print("ChAMP package does not exist, trying to install")
      install.packages("ChAMP")
      if(require("ChAMP")){
        print("ChAMP package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install ChAMP package at: https://www.bioconductor.org/packages/release/bioc/html/ChAMP.html")
      }
    }
    #Load the arraytype and beta matrix
    arraytype = arraytype
    beta=as.matrix(beta)
    ProbeLasso <- champ.DMR(beta = beta,pheno=pd$Sample_Group,method="ProbeLasso", arraytype = arraytype, minProbes = 3)
    ProbeLasso=ProbeLasso[[1]]
    ProbeLasso$methodname="ProbeLasso"
    ProbeLasso
  }



