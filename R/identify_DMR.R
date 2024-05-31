#' Identify the differentially methylated regions using bumphunter method.
#' This function identifies the differentially methylated regions by using bumphunter method. Users should provide corresponding
#' file and parameters, including the methylation beta matrix file, the phenomenon file with sample information, the type of platform, and the threshold for minimum number of probes in a DMR.
#' @param beta The methylation beta values matrix of samples from two groups.
#' @param method The input DMR identification methods, including bumphunter, comb-p, ipDMR, mCSEA, ProbeLasso, and seqlm. One or more methods are allowed to be used.
#' @param pd The phenomenon data of sample information (With Sample group column named Sample_Group).
#' @param arraytype The type of methylation array, including methylation 450K and EPIC array.
#' @param group1 The group name of group one.
#' @param group2 The group name of group two.
#' @param minProbes The threshold for minimum number of probes in a DMR.
#' @param regionsTypes The predefined regions to be analyzed. NULL to skip this step and use customAnnotation (This parameter is only required by mCSEA).
#' @return The DMR set detected by input DMR identification methods.
#' @export
#' @examples
#' ## For 450K data:
#' beta = readRDS(system.file("extdata", "beta_450K.RDS", package = 'DMRIntTk'))
#' pd = read.csv(system.file("extdata", "pd_450K.csv", package = 'DMRIntTk'))
#' totalDMR = identify_DMR(beta = beta, method = c("bumphunter","combp","ipDMR","mCSEA","ProbeLasso","seqlm"), pheno = pd, arraytype = "450K", group1 = "Tumor", group2 = "Normal", minProbes = 3, regionsTypes = "promoter")
#'
#' ##For EPIC data:
#' beta = readRDS(system.file("extdata", "beta_EPIC.RDS", package = 'DMRIntTk'))
#' pd = read.csv(system.file("extdata", "pd_EPIC.csv", package = 'DMRIntTk'))
#' totalDMR = identify_DMR(beta = beta, method = c("bumphunter","combp","ipDMR","mCSEA","seqlm"), pheno = pd, arraytype = "EPIC", group1 = "CRC", group2 = "Polyp", minProbes = 3)
identify_DMR<-
  function(beta = beta, pheno = pd , method = c("bumphunter","combp","ipDMR","mCSEA","ProbeLasso","seqlm"),  arraytype = c("450K", "EPIC"), group1 = group1, group2 = group2, minProbes = minProbes,
           regionsTypes = c("promoters","genes","CGI")){
    library(BiocManager)
    arraytype = match.arg(arraytype)
    regionsTypes = match.arg(regionsTypes)
    beta = as.matrix(beta)
    library(dplyr)
    if("bumphunter" %in% method | "ProbeLasso" %in% method){
    if(require("ChAMP")){
      print("ChAMP package is successfully loaded.")
      library(ChAMP)
    } else {
      print("ChAMP package does not exist, trying to install")
      BiocManager::install("ChAMP")
      if(require("ChAMP")){
        print("ChAMP package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install ChAMP package at: https://www.bioconductor.org/packages/release/bioc/html/ChAMP.html")
      }
    }
    if("bumphunter" %in% method){
      message("-----------------------------")
      message("[<<<<< BUMPHUNTER START >>>>>]")
      message("-----------------------------")

      bumphunter <- champ.DMR(beta = beta,pheno=pd$Sample_Group, method="Bumphunter", arraytype = arraytype, minProbes = minProbes)
      bumphunter=bumphunter[[1]]
      bumphunter$methodname="bumphunter"
      bumphunter<- bumphunter %>% dplyr::select(seqnames, start, end, methodname)
      colnames(bumphunter)[1]="chr"
      message("-----------------------------")
      message("[<<<<< BUMPHUNTER END >>>>>]")
      message("-----------------------------")
    }
    if("ProbeLasso" %in% method){
      message("-----------------------------")
      message("[<<<<< PROBELASSO START >>>>>]")
      message("-----------------------------")
      ProbeLasso <- champ.DMR(beta = beta,pheno=pd$Sample_Group,method="ProbeLasso", arraytype = arraytype, minProbes = 3)
      ProbeLasso=ProbeLasso[[1]]
      ProbeLasso$methodname="ProbeLasso"
      ProbeLasso<- ProbeLasso %>% dplyr::select(seqnames, start, end, methodname)
      colnames(ProbeLasso)[1]="chr"
      message("-----------------------------")
      message("[<<<<< PROBELASSO END >>>>>]")
      message("-----------------------------")
    }
  }
    if("combp" %in% method | "ipDMR" %in% method){
    if(require("ENmix")){
      print("ENmix package is successfully loaded.")
      library(ENmix)
    } else {
      print("ENmix package does not exist, trying to install")
      BiocManager::install("ENmix")
      if(require("ENmix")){
        print("ENmix package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install ENmix package at: https://www.bioconductor.org/packages/release/bioc/html/ENmix.html")
      }
    }
    if(require("CpGassoc")){
      print("CpGassoc package is successfully loaded.")
      library(CpGassoc)
    } else {
      print("CpGassoc package does not exist, trying to install")
      BiocManager::install("CpGassoc")
      if(require("CpGassoc")){
        print("CpGassoc package is successfully installed, loading now")
      } else {
        stop("Installation failed, please go to Bioconductor to install CpGassoc package at: https://www.bioconductor.org/packages/release/bioc/html/CpGassoc.html")
      }
    }
    if(arraytype == "450K"){
      information <- DMRIntTk0411::data_450K
    }else{
      information <- DMRIntTk0411::data_EPIC
    }
      pheno=pd$Sample_Group
      pheno[pheno==group1]<-0
      pheno[pheno==group2]<-1
      pheno <- as.matrix(pheno)
      results<-cpg.assoc(beta,pheno,large.data=FALSE)
      probcg <- rownames(beta)
      cpg <- information[which(information$Name %in% probcg), ]
      combp_beta <- as.data.frame(matrix(NA, nrow(cpg), 5))
      colnames(combp_beta) <- c("chr", "start", "end", "p", "probe")
      combp_beta$chr <- cpg$CHR
      combp_beta$start <- cpg$MAPINFO
      combp_beta$end <- cpg$MAPINFO + 1
      combp_beta$p <- results$results$P.value
      combp_beta$probe <- cpg$Name
    if("combp" %in% method){
      message("-----------------------------")
      message("[<<<<< COMB-P START >>>>>]")
      message("-----------------------------")
      combp<-combp(combp_beta,seed = 0.05,region_plot = F,mht_plot = F,nCores = 1)
      combp<-read.csv("resu_combp.csv")
      combp<-combp[combp$fdr<0.05&combp$nprobe>minProbes,]
      combp$methodname="combp"
      combp$chr=paste0("chr",combp$chr)
      combp<- combp %>% dplyr::select(chr, start, end, methodname)
      message("-----------------------------")
      message("[<<<<< COMB-P END >>>>>]")
      message("-----------------------------")
    }
    if("ipDMR" %in% method){
      message("-----------------------------")
      message("[<<<<< IPDMR START >>>>>]")
      message("-----------------------------")
      ipDMR<-ipdmr(combp_beta,seed = 0.05,region_plot = F ,mht_plot = F)
      ipDMR<-read.csv("resu_ipdmr.csv")
      ipDMR<-ipDMR[ipDMR$fdr<0.05&ipDMR$nprobe>minProbes,]
      ipDMR$methodname="ipDMR"
      ipDMR$chr=paste0("chr",ipDMR$chr)
      ipDMR<- ipDMR %>% dplyr::select(chr, start, end, methodname)
      message("-----------------------------")
      message("[<<<<< IPDMR END >>>>>]")
      message("-----------------------------")
    }
  }
    if("mCSEA" %in% method){
      if(require("mCSEA")){
        print("mCSEA package is successfully loaded.")
        library(mCSEA)
      } else {
        print("mCSEA package does not exist, trying to install")
        install.packages("mCSEA")
        if(require("mCSEA")){
          print("mCSEA package is successfully installed, loading now")
        } else {
          stop("Installation failed, please go to Bioconductor to install mCSEA package at: https://www.bioconductor.org/packages/release/bioc/html/mCSEA.html")
        }
      }
      message("-----------------------------")
      message("[<<<<< MCSEA START >>>>>]")
      message("-----------------------------")
      if(arraytype == "450K"){
        platform <- "450k"
        information <- DMRIntTk0411::data_450K
      }else{
        platform <- "EPIC"
        information <- DMRIntTk0411::data_EPIC
      }
      Mv1 <- log2(beta/(1-beta))
      pheno <- as.data.frame(pd$Sample_Group)
      cnames <- colnames(beta)
      rownames(pheno) <- cnames
      myRank <- rankProbes(Mv1, pheno, caseGroup = group1 ,refGroup = group2, typeInput = 'M', typeAnalysis = 'M')

      mCSEA <- mCSEATest(myRank, beta, pheno,
                         regionsTypes = regionsTypes, platform = platform, minCpGs = minProbes)
      mCSEA=mCSEA[[regionsTypes]]
      mCSEA$NES=as.numeric(mCSEA$NES)
      mCSEA<-mCSEA[abs(mCSEA$NES)>1 & mCSEA$padj<0.05,]
      mCSEA$chr=""
      mCSEA$start=""
      mCSEA$end=""
      mCSEA$probe_length=""
      for (i in 1:dim(mCSEA)[1]) {
        probes=strsplit(mCSEA$leadingEdge[i],split=", ")[[1]]
        mCSEA$chr[i]=paste0("chr",information[which(information$Name==probes[1]),]$CHR)
        mCSEA$start[i]=information[which(information$Name==probes[1]),]$MAPINFO
        mCSEA$end[i]=information[which(information$Name==probes[1]),]$MAPINFO
        mCSEA$probe_length[i]=length(probes)
        if(mCSEA$probe_length[i]>1){
        for (j in 2:length(probes)) {
        tempstart=information[which(information$Name==probes[j]),]$MAPINFO
        tempend=information[which(information$Name==probes[j]),]$MAPINFO
        if(tempstart<mCSEA$start[i]){
          mCSEA$start[i]=tempstart
        }
        if(tempend>mCSEA$end[i]){
          mCSEA$end[i]=tempend
          }
         }
        }
      }
      mCSEA=mCSEA[mCSEA$probe_length>=minProbes,]
      mCSEA$methodname="mCSEA"
      mCSEA<- mCSEA %>% dplyr::select(chr, start, end, methodname)
      message("-----------------------------")
      message("[<<<<< MCSEA END >>>>>]")
      message("-----------------------------")
    }
    if("seqlm" %in% method){
      if(require("devtools")){
        print("devtools package is successfully loaded.")
        library(devtools)
      } else {
        print("devtools package does not exist, trying to install")
        BiocManager::install("devtools")
        if(require("devtools")){
          print("devtools package is successfully installed, loading now")
        } else {
          stop("Installation failed, please go to Bioconductor to install devtools package at: https://www.bioconductor.org/packages/release/bioc/html/devtools.html")
        }
      }
      if(require("seqlm")){
        print("seqlm package is successfully loaded.")
        library(seqlm)
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
        BiocManager::install("GenomicRanges")
        if(require("GenomicRanges")){
          print("GenomicRanges package is successfully installed, loading now")
        } else {
          stop("Installation failed, please go to Bioconductor to install GenomicRanges package at: https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html")
        }
      }
      message("-----------------------------")
      message("[<<<<< SEQLM START >>>>>]")
      message("-----------------------------")
    if(arraytype == "450K"){
      Gr1<-DMRIntTk0411::Gr1_450K
    }else{
      Gr1<-DMRIntTk0411::Gr1_EPIC
    }
     lsegments <-  seqlm(values = beta,
                        genome_information = Gr1, annotation =  factor(pd$Sample_Group))
     sig<-lsegments[lsegments$fdr<0.05]
     seqlm=sig[sig$length>1,]
     write.csv(seqlm,"seqlm.csv",row.names = F)
     seqlm<-read.csv("seqlm.csv")
     seqlm$nprobe=""
     seqlm_probe=apply(seqlm ,1,function(data,probe){
      n= length(strsplit(data[probe],";")[[1]])
      return(n)
    },probe='probes')
     seqlm$nprobe=seqlm_probe
     seqlm<-seqlm[seqlm$nprobe>minProbes,]
     seqlm$methodname="seqlm"
     seqlm$seqnames=paste0("chr",seqlm$seqnames)
     seqlm<- seqlm %>% dplyr::select(seqnames, start, end, methodname)
     colnames(seqlm)[1]="chr"
     message("-----------------------------")
     message("[<<<<< SEQLM END >>>>>]")
     message("-----------------------------")

    }
  totalDMR=data.frame()
  library(plyr)
  for(i in 1:length(method)){
  temp=get(method[i])
  totalDMR=rbind.fill(totalDMR,temp)
  colnames(totalDMR)[1]="chr"
  }
  totalDMR=totalDMR[totalDMR$probe_length>minProbes,]
  totalDMR
}

