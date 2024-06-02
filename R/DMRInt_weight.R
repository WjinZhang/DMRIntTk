#' Calculate the weight of intervals.
#'
#' @param bin The pre-split genomic intervals with calculated covered DMRs (obtained from DMRInt_method function).
#' @param input_DMR The DMR set (obtained from the function DMRInt_input function).
#' @param pd The phenomenon data of sample information (With Sample group column named Sample_Group).
#' @param beta  The methylation beta values matrix of samples from two groups.
#' @param group1 The group name of group one.
#' @param group2 The group name of group two.
#' @return The bins with calculated weights.
#' @export
#'
#' @examples
#' bin_weight=DMRInt_weight(bin_method, input_DMR, pd, beta, group1, group2)
DMRInt_weight=function(bin, input_DMR, pd, beta, group1, group2){
  t = max(abs(input_DMR$diff))
  count = 0
  while(t < 1){
    t = t*10
    count = count + 1
}
    calculate_DMRn <- function(input_DMR,count){
      message(sprintf("[%s] # Calculate the DMR score at all thresholds",
                      Sys.time()))
      input_DMR$start=as.numeric(input_DMR$start)
      input_DMR$end=as.numeric(input_DMR$end)
      interval <- as.data.frame(matrix(NA, 10, 2))
      rownames(interval) <- seq(0.1^count,10*0.1^count,0.1^count)
      colnames(interval) <- c("Tprb","Tlength")
      interval$Tprb <- 0
      interval$Tlength <- 0
      message("Calculate the total number of cpg, probe and DMR length for each interval...")
      for(n in c(1:10)){
        t1 <- (n-1)*(0.1^count)
        t2 <- n*(0.1^count)
        tmp <- input_DMR[which(abs(input_DMR$diff) > t1 & abs(input_DMR$diff) <= t2), ]
        interval[n, 1] <- sum(tmp$probe_length)
        interval[n, 2] <- sum(tmp[, "end"]-tmp[, "start"])

      }
      w <- seq(0.1^count/2, (9*0.1^count+10*0.1^count)/2, 0.1^count)
      res_Qn <- c()
      for (i in c(1:10)) {
        interval$weight <- seq(0.1^count/2, (9*0.1^count+10*0.1^count)/2, 0.1^count)
        tmp <- interval[i:dim(interval)[1], ]
        tmp <- tmp[tmp$Tprb != 0, ]
        if (dim(tmp)[1] == 0) {
          res_Qn <- append(res_Qn, 0)
          next
        }
        Qn_tmp <- sum(tmp$weight * tmp$Tprb^4/tmp$Tlength)/sum(w[i:10])
        Qn_tmp <- log2(Qn_tmp + 1)
        res_Qn <- append(res_Qn, Qn_tmp)
      }
      message(sprintf("[%s] Done", Sys.time()))
      res_Qn
    }

    method=unique(input_DMR$methodname)
    #Construct the score matrix
    DMRn_matrix=as.data.frame(matrix(data=NA,ncol = 10,nrow = length(method)))
    colnames(DMRn_matrix)=paste0("t=",as.character(seq(0,1*0.1^(count-1),0.1^count)[1:10]))

    row.names(DMRn_matrix)=method

    for(name in method){
      temp=input_DMR[which(input_DMR$methodname==name),]
      DMRn=calculate_DMRn(temp,count)
      DMRn_matrix[which(row.names(DMRn_matrix)==name),]=DMRn
    }
    weight_m<-as.data.frame(matrix(data = NA,ncol = 10,nrow = length(method)))
    colnames(weight_m)<-colnames(DMRn_matrix)
    row.names(weight_m)<-row.names(DMRn_matrix)
    for(i in 1:length(method)){
      weight=c()
      for(j in 1:10){
        tempw=DMRn_matrix[i,j]/sum(DMRn_matrix[,j])
        weight=append(weight,tempw)
      }
      weight_m[i,]=weight
    }
    cat(sprintf("[%s] Done\n", Sys.time()))

probe_diff=apply(beta[,which(pd$Sample_Group==group1)],1,mean)-apply(beta[,which(pd$Sample_Group==group2)],1,mean)
cg = strsplit(bin$probe,', ')
bin_diff=sapply(1:nrow(bin), function(i){
  if(i %% 1000 == 0){print(i)}
  combine_cg = intersect(cg[[i]],rownames(beta))
  return(mean(probe_diff[combine_cg]))
})

bin$diff=bin_diff
bin=bin[!is.na(bin$diff),]

#Calculate weight for each bin
bin$weight_count=""
t = max(abs(bin$diff))
count = 0
while(t < 1){
  t = t*10
  count = count + 1
}
#maxdiff defines the last column with values in the weight matrix.
maxdiff=as.numeric(gsub("t=","",colnames(weight_m[which(is.na(weight_m[1,]))])[1]))-0.1^count
#Calculate the sum of weights of all methods covering the bin under the corresponding diff
for(i in 1:dim(bin)[1]){
  tempdiff=round((abs(bin$diff[i])-0.5*(0.1^count)),count)
  if(bin$method_count[i]==""){
    bin$weight_count[i]=0
  }else if(tempdiff<maxdiff){
    method=unique(strsplit(bin$method[i],split = ";")[[1]])
    weight_count=0
    for(j in 1:length(method)){
      temp=weight_m[which(row.names(weight_m)==method[j]),(tempdiff+0.1^count)/0.1^count]
      weight_count=weight_count+temp
    }
    bin$weight_count[i]=weight_count
  }else{
    method=unique(strsplit(bin$method[i],split = ";")[[1]])
    weight_count=0
    for(j in 1:length(method)){
      temp=weight_m[which(row.names(weight_m)==method[j]),which(is.na(weight_m[1,]))[1]-1]
      weight_count=weight_count+temp
    }
    bin$weight_count[i]=weight_count
  }
}

bin$weight_count=as.numeric(bin$weight_count)
bin$weight_diff=bin$weight_count*abs(bin$diff)
    bin
}
