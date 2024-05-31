#' Prepare the input of DMRIntTk.
#' This function finds the probes included in each DMR, and calculates the methylation differences of probes and DMRs.
#' @param totalDMR Multiple DMR sets predicted from different methods.
#' @param beta  The methylation beta values matrix of samples from two groups.
#' @param group1 The group name of group one.
#' @param group2 The group name of group two.
#' @param arraytype The type of methylation array, including 450K and EPIC array.
#' @return DMR sets matrix with calculated methylation differences and probes.
#' @export
#' @examples
#' For 450K data:
#' totalDMR = DMRInt_input(totalDMR, beta , group1 = "Tumor", group2 = "Normal" , arraytype = "450K")
#'
#' For EPIC data:
#' totalDMR = DMRInt_input(totalDMR, beta , group1 = "CRC", group2 = "Polyp" , arraytype = "EPIC")

DMRInt_input<-
         function(totalDMR, beta, group1 = group1, group2 = group2, arraytype=c("450K","EPIC")){
         totalDMR$chr=gsub("chr","",totalDMR$chr)
         totalDMR$chr = as.character(totalDMR$chr)
         arraytype=match.arg(arraytype)
         if(arraytype=="450K"){
           illumina=DMRIntTk0411::data_450K
         }else if(arraytype=="EPIC"){
           illumina=DMRIntTk0411::data_EPIC
         }else{
           stop("The arraytype must be illumina 450K or EPIC...")
         }

         for(i in 1:dim(totalDMR)[1]){
           temp=illumina[illumina$CHR==as.character(totalDMR[i,]$chr),]
           index <- temp[as.integer(totalDMR[i,]$start) <= temp$MAPINFO & as.integer(totalDMR[i,]$end) >=
                          temp$MAPINFO, ]
           probe <- paste(index$Name, collapse = ", ")
           length <- length(index$Name)

           totalDMR$probe[i]=probe
           totalDMR$probe_length[i]=length
         }
         beta=as.data.frame(beta)
         probe_diff=apply(beta[,which(pd$Sample_Group==group1)],1,mean)-apply(beta[,which(pd$Sample_Group==group2)],1,mean)
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



