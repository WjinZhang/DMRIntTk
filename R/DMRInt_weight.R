#' Calculate the weight of intervals.
#'
#' @param interval The pre-split genomic intervals with calculated covered DMRs(obtained from DMRIntTk_method function).
#' @param weight_m The weight matrix obtained from DMRInt_matrix function.
#'
#' @return The interval with calculated weights.
#' @export
#'
#' @examples
#' interval_weight=DMRInt_weight(interval_method,weight_m,beta,case,control)
DMRInt_weight=function(interval,weight_m,beta,case,control){

#Calculate the methylation difference of each interval based on the methylation difference of probes
#than included in the interval
probe_diff=apply(beta[,case],1,mean)-apply(beta[,control],1,mean)
cg = strsplit(interval$probe,', ')
interval_diff=sapply(1:dim(interval)[1], function(i){
  if(i %% 1000 == 0){print(i)}
  combine_cg = intersect(cg[[i]],rownames(beta))
  return(mean(probe_diff[combine_cg]))
})
interval$diff=interval_diff
interval=interval[!is.na(interval$diff),]

#Calculate weight for each interval based on its methylation difference and the sum of reliability
#of all methods that cover the interval
interval$weight_count=""
t = max(abs(interval$diff))
count = 0
while(t < 1){
  t = t*10
  count = count + 1
}
#maxdiff defines the last column with values(not NA) in the weight matrix.
maxdiff=as.numeric(gsub("t=","",colnames(weight_m[which(is.na(weight_m[1,]))])[1]))-0.1^count
#Calculate the sum of weights of all methods covering the interval under the corresponding methylation difference
for(i in 1:dim(interval)[1]){
  tempdiff=round((abs(interval$diff[i])-0.5*(0.1^count)),count)
  if(interval$methodcount[i]==0){
    interval$weight_count[i]=0
  }else if(tempdiff<maxdiff){
    method=unique(strsplit(interval$method[i],split = ";")[[1]])
    weight_count=0
    for(j in 1:length(method)){
      temp=weight_m[which(row.names(weight_m)==method[j]),(tempdiff+0.1^count)/0.1^count]
      weight_count=weight_count+temp
    }
    interval$weight_count[i]=weight_count
  }else{
    method=unique(strsplit(interval$method[i],split = ";")[[1]])
    weight_count=0
    for(j in 1:length(method)){
      temp=weight_m[which(row.names(weight_m)==method[j]),which(is.na(weight_m[1,]))[1]-1]
      weight_count=weight_count+temp
    }
    interval$weight_count[i]=weight_count
  }
}
#Multiply the reliability weights of the interval by the its methylation difference
#as the final weight of the interval
interval$weight_count=as.numeric(interval$weight_count)
interval$weight_diff=interval$weight_count*abs(interval$diff)
    interval
}
