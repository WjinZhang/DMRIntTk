#' Construct the reliability weight matrix.
#' This function calculates the DMRscore and weights of all methods under different methylation difference thresholds and constructs the reliability weight matrix for each method.
#' @param totalDMR DMRs sets obtained from the function preinput.
#' @return The reliability weight matrix.
#' @export
#'
#' @examples
#' weight_m = DMRInt_matrix(totalDMR)
#
DMRInt_matrix <-
   function(totalDMR){
   # Calculate DMRn of each method under different thresholds
   t = max(abs(totalDMR$diff))
   count = 0
   while(t < 1){
   t = t*10
  count = count + 1
  }

  calculate_DMRn <- function(totalDMR,count){
    message(sprintf("[%s] # Calculate the DMR score at all thresholds",
                    Sys.time()))
    interval <- as.data.frame(matrix(NA, 10, 2))
    rownames(interval) <- seq(0.1^count,10*0.1^count,0.1^count)
    colnames(interval) <- c("Tprb","Tlength")
    interval$Tprb <- 0
    interval$Tlength <- 0
    message("Calculate the total number of cpg, probe and DMR length for each interval...")
    for(n in c(1:10)){
      t1 <- (n-1)*(0.1^count)
      t2 <- n*(0.1^count)
      tmp <- totalDMR[which(abs(totalDMR[, 7]) > t1 & abs(totalDMR[, 7]) <= t2), ]
      interval[n, 1] <- sum(tmp[, 6])
      interval[n, 2] <- sum(tmp[, 3]-tmp[, 2])

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
  method=unique(totalDMR$methodname)
  #Construct the DMRn matrix
  DMRn_matrix=as.data.frame(matrix(data=NA,ncol = 10,nrow = length(method)))
  colnames(DMRn_matrix)=paste0("t=",as.character(seq(0,1*0.1^(count-1),0.1^count)[1:10]))

  row.names(DMRn_matrix)=method

  for(name in method){
    temp=totalDMR[which(totalDMR$methodname==name),]
    DMRn=calculate_DMRn(temp,count)
    DMRn_matrix[which(row.names(DMRn_matrix)==name),]=DMRn
  }
  #Construct the weight matrix by normalization
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
   weight_m
}




