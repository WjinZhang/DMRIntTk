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
#' totalDMR = read.csv(system.file("extdata","totalDMR.csv",package = 'DMRIntTk'))
#' beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
#' case = 1:5
#' control = 6:10
#' totalDMR=DMRInt_input(totalDMR, beta, case,control, arraytype = "450K")
#'
     DMRInt_input<- 
               function(totalDMR, beta, case_pos=case, con_pos=control, arraytype=c("450K","EPIC")){
               totalDMR$chr=gsub("chr","",totalDMR$chr)
               totalDMR$chr = as.character(totalDMR$chr)
               illumina450k <- DMRIntTk::illumina450k_hg19
               illuminaEPIC <- DMRIntTk::illuminaEPIC_hg19
               arraytype=match.arg(arraytype)
               if(arraytype=="450K"){
                 illumina=illumina450k
               }else if(arraytype=="EPIC"){
                 illumina=illuminaEPIC
               }else{
                 stop("The arraytype must be illumina 450K or EPIC...")
               }

               for(i in 1:dim(totalDMR)[1]){
                 temp=illumina[illumina$CHR==as.character(totalDMR[i,]$chr),]
                 index <- temp[as.integer(totalDMR[i,]$start) <= temp$MAPINFO & as.integer(totalDMR[i,]$end) >= 
                                temp$MAPINFO, ]
                 probe <- paste(rownames(index), collapse = ", ")
                 length <- length(rownames(index))

                 totalDMR$probe[i]=probe
                 totalDMR$probe_length[i]=length
               }
               beta=as.data.frame(beta)
               probe_diff=apply(beta[,case],1,mean)-apply(beta[,control],1,mean)
               beta$probe_diff=probe_diff
               cg = strsplit(totalDMR$probe,', ')
               totalDMR_diff=sapply(1:dim(totalDMR)[1], function(i){
                 if(i %% 1000 == 0){print(i)}
                 combine_cg = intersect(cg[[i]],rownames(beta))
                 return(mean(probe_diff[combine_cg]))
               })
               totalDMR$diff=totalDMR_diff
               cat(sprintf("[%s] Done\n", Sys.time()))
               totalDMR
}



