
# Introduction
 **DMRIntTk** is a toolkit for **integrating DMR sets** predicted by different methods on a same methylation array dataset based on density peak clustering algorithm.
 In DMR integration, it contains five main functions including **DMRInt_input**, **DMRInt_matrix**, **DMRInt_method**, **DMRInt_weight** and **DMRInt_densitypeak**.
 Before DMR integration, DMRIntTk provides several state-of-art DMR detection functions that help users obtain DMR sets with ease, including **DMRInt_bumphunter**, **DMRInt_ProbeLasso**, **DMRInt_combp**, **DMRInt_ipDMR**, **DMRInt_mCSEA** and **DMRInt_seqlm**.
 The main functions are as the following picture.

 ![image](https://github.com/WjinZhang/DMRIntTk/blob/main/Workflow.jpg)

A schematic diagram of DMRIntTk. (a) Data pre-processing and DMR identication steps output DMR sets are used as standard input in DMRIntTk. There are four major parts in DMRIntTk, including (b) segmenting the genome, (c) constructing the reliability matrix, (d) weighting bins and (e) integrating bins based on density peak clustering. The example of each part are shown.
 # Getting started 
 These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
 ## Prerequisites
 R (version >= 2.10) is required to be installed before using DMRIntTk.

 ## Installing
 You can input the following commands to install the DMRIntTk R package.
 
 ```R
 install.packages("dplyr")
 install.packages("devtools")
 install.packages("ChAMP")
 install.packages("ENmix")
 install.packages("CpGassoc")
 install.packages("GenomicRanges")
 install.packages("mCSEA")
 install_github("raivokolde/seqlm")
 install_github("WjinZhang/DMRIntTk")
 library(DMRIntTk)
 ```
 ## Running the tests

 ### DMR sets identification
 Since DMR integration requires multiple DMR sets predicted by different methods as inputs, users should either have the self-identified multiple DMR sets, or directly use the DMR detection functions provided by DMRIntTk to identify DMR sets.
 For the latter situation, with the **methylation beta value matrix "beta"** and **the phenomenon information "pd"**, users can easily obtain the desired DMR sets with following functions(p.s.: the arraytype and minimum 
 probes can be customized, here we took 450K array and 3 probes for the example):
 ```R
 beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
 pd = read.csv(system.file("extdata","pd.csv",package = 'DMRIntTk'))
 ```
 #### DMRInt_bumphunter
 ```R
 bumphunter = DMRInt_bumphunter(beta = beta, pheno=pd$Sample_Group, arraytype = "450K", minProbes = 3)
 ```
 #### DMRInt_ProbeLasso
 ```R
 ProbeLasso = DMRInt_ProbeLasso(beta = beta, pheno=pd$Sample_Group, arraytype = "450K", minProbes = 3)
```
 #### DMRInt_combp
 ```R
 combp = DMRInt_combp(beta = beta, pd = pd, arraytype = "450K", minProbes = 3)
```
 #### DMRInt_ipDMR
 ```R
 ipDMR = DMRInt_ipDMR(beta = beta, pd = pd, arraytype = "450K", minProbes = 3)
```
 #### DMRInt_mCSEA 
 ```R
 mCSEA = DMRInt_mCSEA(beta = beta, pd = pd, caseGroup= "Tumor", refGroup = "Normal", regionsTypes = "promoters", platform = "450k", minCpGs = 3)
```
 #### DMRInt_seqlm
 ```R
 seqlm = DMRInt_seqlm(beta = beta, pd = pd, arraytype = "450K", minCpGs = 3)
 ```

 ### DMR sets integration
 For self-identified DMR sets, users should organize them into the total DMR file containing **chromosome, start, end and methodnames**. The example total DMR file is provided.
 ```R
 beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
 totalDMR = read.csv(system.file("extdata","totalDMR.csv",package = 'DMRIntTk'))
 ```
 For the DMR sets identified by DMRIntTk, users can get the total DMR set using following codes:
 ```R
 DMRstring = c("chr","start","end", "methodname")
 totalDMR = rbind(bumphunter[,DMRstring], ProbeLasso[,DMRstring], mCSEA[,DMRstring], seqlm[,DMRstring], combp[,DMRstring], ipDMR[,DMRstring])
 ```
 #### DMRInt_input
 Prepare the input of DMRIntTk.
 This function **finds the probes included in each DMR, and calculates the methylation differences of each probe and DMR**.
 DMRInt_input function needs two files as input : 1. **Total DMR sets** obtained from different methods and 2. **Methylation level beta matrix** of all samples. Plus, "case" and "control" are the columm numbers of two groups of samples, respectively.
```R
case = grep("Tumor",colnames(beta))
control = grep("Normal",colnames(beta))
totalDMR = DMRInt_input(totalDMR, beta, case, control, arraytype = "450K")
```
 
 #### DMRInt_matrix
 This function **calculates the DMRscore and weights of all methods** under different methylation difference thresholds and
 **constructs the reliability matrix** for each method. It requires the DMR sets file processed by the function DMRInt_input.

 ```R
 weight_m = DMRInt_matrix(totalDMR)
 ```
 
 #### DMRInt_method
This function determines DMRs that cover each bin, and calculate the numbers of methods that cover the bin. 
It requires the **pre-split genomic bins file** and **total DMR sets** processed by the function DMRInt_input.

```R
interval_method = DMRInt_method(totalDMR, arraytype = "450K" )
```
 #### DMRInt_weight
 This function **calculates the weights of bins**. It requires **the pre-split genomic bins file**(processed from DMRIntTk_method function) 
 and **the reliability matrix** obtained from DMRInt_matrix function.
 
 
```R
interval_weight = DMRInt_weight(interval_method,weight_m,beta,case,control)
```
 
 #### DMRInt_densitypeak
 This function **clusters all bins** based on density peak algorithm. It needs the **bins file** with calculated weights(obtained from DMRInt_weight function)
 and **total DMR sets**(obtained from DMRIntTk_input function).

```R
Res = DMRInt_densitypeak(interval_weight, totalDMR, prefer = "probe", arraytype = "450K")
```
 
 # Built With
  R is a free software environment for statistical computing and graphics.
  
 # Authors
* Xiaoqing Peng - [Central South University](https://life.csu.edu.cn/jsxx.jsp?urltype=news.NewsContentUrl&wbtreeid=1815&wbnewsid=3625)
* Wenjin Zhang - [Central South University](https://life.csu.edu.cn/)

# Version update
* **DMRIntTk** v0.1.0 -  News version 0.1.0, 2023.09.07. The first version.
