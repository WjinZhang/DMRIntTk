#' Identify the differentially methylated regions using seqlm method.
#' This function identifies the differentially methylated regions by using seqlm method. Users should provide corresponding
#' file and parameters, including the methylation beta matrix file, the phenomenon file with sample information, the type of platform, and the threshold for minimum number of probes in a DMR.
#' @param beta The methylation beta values matrix of samples from two groups.
#' @param pd The phenomenon data with sample information (With Sample group column named Sample_Group).
#' @param arraytype The type of methylation array, 450K is allowed.
#' @param minProbes The threshold for minimum number of probes in a DMR.
#' @return The DMR set detected by seqlm.
#' @export
#' @examples
#' beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
#' pd=read.csv(system.file("extdata","pd.csv",package = 'DMRIntTk'))
#' seqlm=DMRInt_seqlm(beta = beta, pd = pd, arraytype = "450K", minCpGs = 3)
#'
DMRInt_seqlm<-
  function(beta, pd, arraytype="450K", minCpGs){
    #Check the installation of seqlm and related packages(seqlm can only be installed by github tools "devtools")

    if(require("devtools")){
      print("devtools package is successfully loaded.")
    } else {
      print("devtools package does not exist, trying to install")
      install.packages("devtools")
      if(require("devtools")){
        print("devtools package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install devtools package at: https://www.bioconductor.org/packages/release/bioc/html/devtools.html")
      }
    }
    if(require("seqlm")){
      print("seqlm package is successfully loaded.")
    } else {
      print("seqlm package does not exist, trying to install")
      library(devtools)
      install_github("raivokolde/seqlm")
      if(require("seqlm")){
        print("seqlm package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install seqlm package at: https://github.com/raivokolde/seqlm")
      }
    }
    if(require("GenomicRanges")){
      print("GenomicRanges package is successfully loaded.")
    } else {
      print("GenomicRanges package does not exist, trying to install")
      install.packages("GenomicRanges")
      if(require("GenomicRanges")){
        print("GenomicRanges package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install GenomicRanges package at: https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html")
      }
    }
    library(seqlm)
    library(GenomicRanges)
    beta=as.matrix(beta)
    #genome information Gr1 is provided in DMRIntTk package
    Gr1<-DMRIntTk::Gr1
    lsegments <-  seqlm(values = beta,
                        genome_information = Gr1, annotation =  factor(pd$Sample_Group))
    #Remain the significant DMRs with adjusted p value < 0.05
    sig<-lsegments[lsegments$fdr<0.05]
    seqlm=sig[sig$length>1,]
    write.csv(seqlm,"seqlm.csv",row.names = F)
    seqlm<-read.csv("seqlm.csv")
    #Filter out DMRs that contain less than the specified number of probes.
    seqlm$nprobe=""
    seqlm_probe=apply(seqlm ,1,function(data,probe){
      n= length(strsplit(data[probe],";")[[1]])
      return(n)
    },probe='probes')
    seqlm$nprobe=seqlm_probe
    seqlm<-seqlm[seqlm$nprobe>minCpGs,]
    seqlm$methodname="seqlm"
    seqlm
  }



