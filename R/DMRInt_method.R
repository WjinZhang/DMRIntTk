#' Determine DMRs that cover each bin, and calculate the numbers of methods that cover the bin.
#' @importFrom data.table ":="
#' @param input_DMR DMRs sets obtained from the function DMRInt_input.
#' @param arraytype The type of methylation array, including 450K and EPIC array.
#' @return The bin with covered DMRs and methods.
#' @export
#' @examples
#' bin_method = DMRInt_method(input_DMR, arraytype = "450K" )
DMRInt_method=function(input_DMR, arraytype = c("450K", "EPIC")){
    library(dplyr)
    input_DMR<-arrange(input_DMR,chr,start)
    input_DMR$ID<-paste0("DMR",rownames(input_DMR))
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

    input_DMR=data.table(input_DMR)
    bin=data.table(bin)
    bin[, DMR := ""]
    bin[, method := ""]
    bin[, method_count := ""]

    for (i in 1:dim(bin)[1]) {
     overlap <- input_DMR[chr == bin$chr[i] &
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
