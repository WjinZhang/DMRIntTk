#' Cluster all intervals based on density peak algorithm.
#'
#' @param interval The interval with calculated weights(obtained from DMRInt_weight function).
#' @param totalDMR The DMR sets obtained from DMRIntTk_input function.
#'
#' @return The integrated DMR set.
#' @export 
#'
#' @examples
#' Res=DMRInt_densitypeak(interval_weight, totalDMR, prefer = "probe", arraytype = "450K")
DMRInt_densitypeak=function(interval,totalDMR,prefer=c("probe","length"), arraytype=c("450K", "EPIC")){
library(dplyr)
DMR<-select(interval,c(chr,start,end,weight_diff,diff))
DMR$weight_diff=as.numeric(DMR$weight_diff)
#Set the clustering threshold, first divide the maxdiff into tenths,
#set the center threshold from the 2/10th to the 10th/10th fraction,
#and set the neighbor threshold from the 1/10th to the 9th/10th fraction.
totweight_par=seq(from=max(DMR$weight_diff)*0.1,to=max(DMR$weight_diff),by=max(DMR$weight_diff)*0.1)
cenweight_par=totweight_par[2:10]
neiweight_par=totweight_par[1:9]

#Construct DMR proportion matrix (calculate the proportion of DMRs exceeding 1/3 and 1/2 of diff, used to select the clustering threshold later)
DMR_pro=as.data.frame(matrix(data = NA,nrow = 10,ncol = 2))
colnames(DMR_pro)=c("sum_1","sum_2")
rownames(DMR_pro)=c("2_1","3_1","3_2","4_1","4_2","4_3","5_1","5_2","5_3","5_4")
DMR_stat=DMR_pro
colnames(DMR_stat)=c("length","probe")
maxdiff=max(abs(totalDMR$diff))
threshold1=0.33*maxdiff
threshold2=0.5*maxdiff
#Save the all the clustering results under each threshold in the output list
output=list()
#Density peak
for(th_center in cenweight_par[1:4]){
  for(threshold in neiweight_par[1:4]){
    if(th_center>threshold){
      temp <- DMR
      temp <- temp[order(temp$chr,temp$start),]
      temp$index <- as.numeric(rownames(temp))
      temp <- temp[order(temp$index),]
      local_max <- apply(temp,1,FUN = function(data,chr,start,weight_diff){
        nchr <- as.numeric(data[chr])
        nstart <- as.numeric(data[start])
        nweight_diff <- as.numeric(data[weight_diff])
        if(nweight_diff >= th_center){
          localmax <- -1
        }else{
          index <- which(temp$chr == nchr & temp$weight_diff > nweight_diff)
          if(length(index) > 0){
            localmax <- index[which.min(abs(nstart-temp$start[index]))]
          }else{
            localmax <- -1
          }
        }
        localmax
      },chr = 'chr',start = 'start',weight_diff = 'weight_diff')
      
      temp$local_max <- local_max
      
      new_locmax <- apply(temp,1,FUN = function(data,local_max,index){
        local_max <- as.numeric(data[local_max])
        if(local_max == -1){
          ind <- -1
        }else{
          ind <- as.numeric(data[index])
          while(local_max != -1){
            ind <- local_max
            local_max <- temp$local_max[local_max]
          }
        }
        ind
      },local_max = 'local_max',index = 'index')
      temp$local_max <- new_locmax
      
      fake_cen<-which(temp$local_max==-1 & temp$weight_diff<th_center)
      temp$local_max[fake_cen]=0
      
      
      combined_dmr2 <- as.data.frame(matrix(NA,0,7))
      cen_i <- which(temp$local_max == -1)
      for(i in 1:length(cen_i)){
        cen <- as.numeric(cen_i[i])
        nchr <- temp$chr[cen]
        left_boundary <- cen
        while (temp$chr[left_boundary] == nchr & temp$weight_diff[left_boundary] >=threshold
               & ((temp$local_max[left_boundary] == cen) | (temp$local_max[left_boundary] == -1))) {
          if(temp$local_max[left_boundary] == -1) {cen <- temp$index[left_boundary]}
          left_boundary <- left_boundary-1
        }
        cen <- as.numeric(cen_i[i])
        right_boundary <- cen
        
        while(temp$chr[right_boundary] == nchr & temp$weight_diff[right_boundary] >=threshold &
              ((temp$local_max[right_boundary] == cen) | (temp$local_max[right_boundary] == -1))){
          if((temp$local_max[right_boundary] == -1)){cen <- temp$index[right_boundary]}
          right_boundary <- right_boundary+1
          if(right_boundary==(dim(temp)[1]+1)){
            break;
          }
          
        }
        start <- temp$start[left_boundary+1]
        end <- temp$end[right_boundary-1]
        freq <- mean(temp$weight_diff[(left_boundary+1):(right_boundary-1)])
        dif <- mean(temp$diff[(left_boundary+1):(right_boundary-1)])
        
        combined_dmr2 <- rbind(combined_dmr2,c(nchr,start,end,freq,dif,left_boundary+1,right_boundary-1))
      }
      colnames(combined_dmr2) <- c('chr','start','end','weight_diff','diff','left','right')
      combined_dmr2 <- unique(combined_dmr2)
      
      start <- which(diff(combined_dmr2$start) == 0)
      if(length(start) != 0){
        combined_dmr2 <- combined_dmr2[-(start),]}
      
      end <- which(diff(combined_dmr2$end) == 0)
      if(length(end) != 0){
        
        combined_dmr2 <- combined_dmr2[-(end+1),]}
      
      
      if(dim(combined_dmr2)[1]>1){
        colnames(combined_dmr2)[4]<-"weight_diff"
        finaldmr<-as.data.frame(matrix(NA,0,5))
        colnames(finaldmr)<-c("chr","start","end","weight_diff","diff")
        j=1
        for(i in 1:(dim(combined_dmr2)[1]-1)){
          start=combined_dmr2[i,]$start
          chr=combined_dmr2[i,]$chr
          
          #Merge adjacent intervals
          if(combined_dmr2[i,]$end==combined_dmr2[i+1,]$start-1){
            print(i)
            while(combined_dmr2[i,]$chr==combined_dmr2[i+1,]$chr&combined_dmr2[i,]$end==combined_dmr2[i+1,]$start-1 ){
              i=i+1
            }
            end=combined_dmr2[i,]$end
            pos=seq(i:j)
            weight_diff=mean(combined_dmr2[pos,]$weight_diff)
            diff=mean(combined_dmr2[pos,]$diff)
            finaldmr[j,]$chr= chr
            finaldmr[j,]$start=start
            finaldmr[j,]$end=end
            finaldmr[j,]$weight_diff=weight_diff
            finaldmr[j,]$diff=diff
            j=j+1
          }
          else{
            finaldmr[j,]$chr= combined_dmr2[i,]$chr
            finaldmr[j,]$start= combined_dmr2[i,]$start
            finaldmr[j,]$end= combined_dmr2[i,]$end
            finaldmr[j,]$weight_diff= combined_dmr2[i,]$weight_diff
            finaldmr[j,]$diff= combined_dmr2[i,]$diff
            j=j+1
          }
        }
        finaldmr[j,]=combined_dmr2[dim(combined_dmr2)[1],1:5]
        pos=which(!duplicated(finaldmr$end))
        final<-finaldmr[pos,]
        pro1=length(which(abs(final$diff)>threshold1))/dim(final)[1]
        pro2=length(which(abs(final$diff)>threshold2))/dim(final)[1]
        DMR_pro[rownames(DMR_pro)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$sum_1=pro1
        DMR_pro[rownames(DMR_pro)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$sum_2=pro2
        final$chr=as.character(final$chr)
        
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
        
        for(i in 1:dim(final)[1]){
          temp=illumina[illumina$CHR==as.character(final[i,]$chr),]
          index <- temp[as.integer(final[i,]$start) <= temp$MAPINFO & as.integer(final[i,]$end) >= 
                          temp$MAPINFO, ]
          probe <- paste(rownames(index), collapse = ", ")
          length <- length(rownames(index))
          final$probe[i]=probe
          final$probe_length[i]=length
        }
        DMR_stat[rownames(DMR_stat)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$length=sum(final$length)
        DMR_stat[rownames(DMR_stat)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$probe=sum(final$probeLength)
        print(paste0("center=",which(totweight_par==th_center)," neighbor=",which(totweight_par==threshold)))
        output=append(output,list(final))
      }else{  #The case of only 1 interval left after clustering (it may occur when the center threshold is equal to 10/10)
        pro1=length(which(abs(combined_dmr2$diff)>threshold1))/dim(combined_dmr2)[1]
        pro2=length(which(abs(combined_dmr2$diff)>threshold2))/dim(combined_dmr2)[1]
        DMR_pro[rownames(DMR_pro)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$sum_1=pro1
        DMR_pro[rownames(DMR_pro)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$sum_2=pro2
        
        for(i in 1:dim(combined_dmr2)[1]){
          temp=illumina[illumina$CHR==as.character(combined_dmr2[i,]$chr),]
          index <- temp[as.integer(combined_dmr2[i,]$start) <= temp$MAPINFO & as.integer(combined_dmr2[i,]$end) >= 
                          temp$MAPINFO, ]
          probe <- paste(rownames(index), collapse = ", ")
          length <- length(rownames(index))
          combined_dmr2$probe[i]=probe
          combined_dmr2$probe_length[i]=length
        }
        DMR_stat[rownames(DMR_stat)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$length=sum(combined_dmr2$length)
        DMR_stat[rownames(DMR_stat)==paste0(which(totweight_par==th_center),"_",which(totweight_par==threshold)),]$probe=sum(combined_dmr2$probeLength)
        output=append(output,list(combined_dmr2))
      }
    }
  }
}

method=unique(totalDMR$methodname)
sin_pro=as.data.frame(matrix(data = NA,nrow = length(method),ncol = 2))
colnames(sin_pro)=c("sum_1","sum_2")
rownames(sin_pro)=method

for(name in method){
  print(name)
  pro1=length(which(abs(totalDMR[totalDMR$methodname==name,]$diff)>threshold1))/length(which(totalDMR$methodname==name))
  pro2=length(which(abs(totalDMR[totalDMR$methodname==name,]$diff)>threshold2))/length(which(totalDMR$methodname==name))
  sin_pro[which(rownames(sin_pro)==name),]$sum_1=pro1
  sin_pro[which(rownames(sin_pro)==name),]$sum_2=pro2
}

#Determine the rownames of the integrated results whose sum_1 and sum_2 are greater than those of all individual methods
rownames_DMR=rownames(DMR_pro[which(DMR_pro$sum_1 >max(sin_pro$sum_1) & DMR_pro$sum_2>max(sin_pro$sum_2)),])
#Sort all qualified integrated DMR sets by total length and number of probes
length_rank=arrange(DMR_stat[rownames(DMR_stat) %in% rownames_DMR,],desc(length))
probe_rank=arrange(DMR_stat[rownames(DMR_stat) %in% rownames_DMR,],desc(probe))
#Return the recommended DMR set according user's preference(DMR length is the default)
prefer=match.arg(prefer)
if(prefer=="length"){
  Recommend=rownames(length_rank)[1]
}else if(prefer=="probe"){
  Recommend=rownames(probe_rank)[1]
}else{
  stop("The preference must be probe or length...")
}

print(paste0(Recommend," is recommended"))
write.csv(output[[which(rownames(DMR_pro)==Recommend)]],"Integrated DMR set.csv")
output[[which(rownames(DMR_pro)==Recommend)]]
}