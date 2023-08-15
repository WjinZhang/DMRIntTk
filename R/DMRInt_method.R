#' Determine DMRs that cover each interval, and calculate the numbers of methods that cover the interval.
#'
#' @param totalDMR DMRs sets obtained from the function DMRInt_input.
#' @param arraytype The type of methylation array, including 450K and EPIC array.
#'
#' @return The interval with covered DMRs and methods.
#' @export
#'
#' @examples
#' interval_method = DMRInt_method(totalDMR, arraytype = "450K" )
DMRInt_method=function(totalDMR, arraytype = c("450K", "EPIC")){
library(dplyr)
totalDMR<-arrange(totalDMR,chr,start)
totalDMR$ID<-paste0("DMR",rownames(totalDMR))
arraytype=match.arg(arraytype)
if(arraytype=="450K"){
  interval=DMRIntTk::interval450K
}else if(arraytype=="EPIC"){
  interval=DMRIntTk::intervalEPIC
}else{
  stop("The array type must be illumina 450K or EPIC.")
}

interval$DMR<-""

for(i in 1:dim(interval)[1]){
  temp=totalDMR[totalDMR$chr==as.character(interval[i,]$chr),]
  index <- temp[which(interval[i,]$end > temp$start & interval[i,]$start < temp$end), ]
  pos=as.numeric(rownames(index))
  length=length(pos)
  for(j in 1:length){
    interval[i,]$DMR<-paste0(interval[i,]$DMR,totalDMR[pos[j],]$ID, sep=";")
  }
}
interval$DMR<-gsub("NA;;","",interval$DMR)
interval$DMR<-gsub(";$","",interval$DMR)

#Find the number and names of methods that cover each interval
interval$DMRcount<-""
interval$method<-""

DMR_count = apply(interval,1,function(data,DMR){
  length = length(strsplit(data[DMR],";")[[1]])
  return(length)
},DMR='DMR')

interval$DMRcount=DMR_count
for(i in 1:dim(interval)[1]){
  for(j in 1:interval[i,]$DMRcount){
    temp=strsplit(interval[i,]$DMR,split = ";")[[1]][j]
    method=totalDMR[which(totalDMR$ID==temp),]$methodname
    interval[i,]$method<-paste0(interval[i,]$method,";",method)
  }
}
interval$method=gsub(";;","",interval$method)
interval$method=gsub("^;","",interval$method)
interval$methodcount=""

method_count = apply(interval,1,function(data,method){
  length = length(unique(strsplit(data[method],";")[[1]]))
  return(length)
},method='method')
interval$methodcount = method_count
  interval
}