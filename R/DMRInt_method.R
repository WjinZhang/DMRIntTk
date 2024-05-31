#' Determine DMRs that cover each bin, and calculate the numbers of methods that cover the bin.
#' @importFrom data.table ":="
#' @param totalDMR DMRs sets obtained from the function DMRInt_input.
#' @param arraytype The type of methylation array, including 450K and EPIC array.
#' @return The bin with covered DMRs and methods.
#' @export
#' @examples
#' bin_method = DMRInt_method(totalDMR, arraytype = "450K" )
DMRInt_method=function(totalDMR, arraytype = c("450K", "EPIC")){
    library(dplyr)
    totalDMR<-arrange(totalDMR,chr,start)
    totalDMR$ID<-paste0("DMR",rownames(totalDMR))
    arraytype=match.arg(arraytype)
    if(arraytype=="450K"){
      bin=DMRIntTk::bin_450K
    }else if(arraytype=="EPIC"){
      bin=DMRIntTk::bin_EPIC
    }else{
      stop("The array type must be illumina 450K or EPIC.")
    }

    #Find the number and names of methods that cover each bin
    library(data.table)

    totalDMR=data.table(totalDMR)
    bin=data.table(bin)
    bin[, DMR := ""]
    bin[, method := ""]
    bin[, method_count := ""]

    for (i in 1:dim(bin)[1]) {
     overlap <- totalDMR[chr == bin$chr[i] &
                          start < bin$end[i] &
                          end > bin$start[i]]

      if (nrow(overlap) > 0) {
        bin$DMR[i] <- paste(overlap$ID, collapse = ";")
        bin$method[i] <- paste(overlap$method, collapse = ";")
        bin$method_count[i] <- length(unique(overlap$method))
      }
    }
      bin
    }
