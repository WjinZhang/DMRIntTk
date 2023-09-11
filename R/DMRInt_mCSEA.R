#' Identify the differentially methylated regions using mCSEA method.
#' This function identifies the differentially methylated regions by using mCSEA method. Users should provide corresponding
#' file and parameters, including the methylation beta matrix file, the phenomenon file with sample information, the caseGroup
#' and the refGroup, the interested regions types, the type of platform, and the threshold for minimum number of CpGs in a DMR.
#' @param beta The methylation beta values matrix of samples from two groups.
#' @param pd The phenomenon data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==file with sample information (With Sample group column named Sample_Group).
#' @param caseGroup  The name of case group in the pd file.
#' @param refGroup The name of control group in the pd file.
#' @param regionsTypes The region types you are interested in. mCSEA will identify the DMRs in that specific type of regions.
#' @param platform The type of methylation array, including 450K and EPIC array.
#' @param minCpGs The threshold for minimum number of CpGs in a DMR.
#' @return The DMR set detected by mCSEA.
#' @export
#' @examples
#' beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
#' pd=read.csv(system.file("extdata","pd.csv",package = 'DMRIntTk'))
#' mCSEA=DMRInt_mCSEA(beta = beta, pd = pd, caseGroup= "Tumor", refGroup = "Normal", regionsTypes = "promoters", platform = "450k", minCpGs = 3)
#'
DMRInt_mCSEA<-
  function(beta, pd, caseGroup, refGroup, regionsTypes=c("promoters","genes","CGI"), platform=c("450k","EPIC"), minCpGs){
    if(require("mCSEA")){
      print("mCSEA package is successfully loaded.")
    } else {
      print("mCSEA package does not exist, trying to install")
      install.packages("mCSEA")
      if(require("mCSEA")){
        print("mCSEA package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install mCSEA package at: https://www.bioconductor.org/packages/release/bioc/html/mCSEA.html")
      }
    }
    library("mCSEA")
    Mv1 <- log2(beta/(1-beta))
    pheno <- as.data.frame(pd$Sample_Group)
    cnames <- colnames(beta)
    rownames(pheno) <- cnames
    # ??phneo?????ĳ?beta????
    myRank <- rankProbes(Mv1, pheno, caseGroup = caseGroup ,refGroup = refGroup, typeInput = 'M', typeAnalysis = 'M')
    # regionsTypes = "promoters"
    mCSEA <- mCSEATest(myRank, beta, pheno,
                           regionsTypes = regionsTypes, platform = platform, minCpGs = minCpGs)
    mCSEA=mCSEA[[regionsTypes]]
    mCSEA$NES=as.numeric(mCSEA$NES)
    mCSEA<-mCSEA[abs(mCSEA$NES)>1 & mCSEA$padj<0.05,]
    mCSEA$methodname="mCSEA"
    mCSEA
  }



