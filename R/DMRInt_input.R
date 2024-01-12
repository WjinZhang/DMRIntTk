#' Prepare the input of DMRIntTk.
#' This function finds the probes included in each DMR, and calculates the methylation differences of probes and DMRs.
#' @param totalDMR Multiple DMR sets predicted from different methods.
#' @param beta  The methylation beta values matrix of samples from two groups.
#' @param case  The column of case group.
#' @param control The column of control group.
#' @param arraytype The type of methylation array, including 450K and EPIC array.
#' @return DMR sets matrix with calculated methylation differences and probes.
#' @export
#' @examples
#' If you have used the DMR detection functions provided by DMRIntTk and
#' obtained different DMR sets, you can combine them to get the total DMR set:
#' DMRstring = c("chr","start","end", "methodname")
#' totalDMR = rbind(bumphunter[,DMRstring], ProbeLasso[,DMRstring], mCSEA[,DMRstring], seqlm[,DMRstring], combp[,DMRstring], ipDMR[,DMRstring])
#'
#' If you already have the total DMR sets matrix, or want to use our example DMR matrix directly:
#' totalDMR = read.csv(system.file("extdata","totalDMR.csv",package = 'DMRIntTk'))
#' beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
#' case = grep("Tumor",colnames(beta))
#' control = grep("Normal",colnames(beta))
#' totalDMR=DMRInt_input(totalDMR, beta, case, control, arraytype = "450K")
#'
     DMRInt_input<-
               function(totalDMR, beta, case_pos=case, con_pos=control, arraytype=c("450K","EPIC")){
               #Uniform the format of the chromosome column
               totalDMR$chr=gsub("chr","",totalDMR$chr)
               totalDMR$chr = as.character(totalDMR$chr)
               #Load the annotation file of probe position for 450K and 850K
               illumina450k <- DMRIntTk::illumina450k_hg19
               illuminaEPIC <- DMRIntTk::illuminaEPIC_hg19
               #Load the arraytype
               arraytype=match.arg(arraytype)
               if(arraytype=="450K"){
                 illumina=illumina450k
               }else if(arraytype=="EPIC"){
                 illumina=illuminaEPIC
               }else{
                 stop("The arraytype must be illumina 450K or EPIC...")
               }
               #Calculate the probes included in each DMR
               for(i in 1:dim(totalDMR)[1]){
                 temp=illumina[illumina$CHR==as.character(totalDMR[i,]$chr),]
                 index <- temp[as.integer(totalDMR[i,]$start) <= temp$MAPINFO & as.integer(totalDMR[i,]$end) >=
                                temp$MAPINFO, ]
                 probe <- paste(rownames(index), collapse = ", ")
                 length <- length(rownames(index))

                 totalDMR$probe[i]=probe
                 totalDMR$probe_length[i]=length
               }
               #Calculate the methylation difference of each probe
               beta=as.data.frame(beta)
               probe_diff=apply(beta[,case],1,mean)-apply(beta[,control],1,mean)
               beta$probe_diff=probe_diff
               cg = strsplit(totalDMR$probe,', ')
               #Calculate the methylation difference of each DMR
               totalDMR_diff=sapply(1:dim(totalDMR)[1], function(i){
                 if(i %% 1000 == 0){print(i)}
                 combine_cg = intersect(cg[[i]],rownames(beta))
                 return(mean(probe_diff[combine_cg]))
               })
               totalDMR$diff=totalDMR_diff
               cat(sprintf("[%s] Done\n", Sys.time()))
               totalDMR
}



