
# Introduction
 **DMRIntTk** is a toolkit for **integrating DMR sets** predicted by different methods on a same methylation array dataset based on density peak clustering algorithm.
 In DMR integration, it contains five main functions including **identify_DMR**, **DMRInt_input**, **DMRInt_method**, **DMRInt_weight**,  and **DMRInt_densitypeak**.
 Before DMR integration, DMRIntTk provides several state-of-art DMR detection functions that help users obtain DMR sets with ease, including bumphunter, ProbeLasso, combp, ipDMR, mCSEA and seqlm. These methods can be easily conducted by the function **identify_DMR**.
 The main functions are as the following picture.

 ![image](https://github.com/WjinZhang/DMRIntTk/blob/main/Workflow.jpg)

A schematic diagram of DMRIntTk. (a) Data pre-processing and DMR identication steps output DMR sets are used as standard input in DMRIntTk. There are four major parts in DMRIntTk, including (b) segmenting the genome, (c) constructing the reliability matrix, (d) weighting bins and (e) integrating bins based on density peak clustering. The example of each part are shown.
 # Getting started 
 These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
 ## Prerequisites
 R (version >= 4.0) is required to be installed before using DMRIntTk.

 ## Installing
 You can input the following commands to install the DMRIntTk R package.
 
 ```R
 library(BiocManager)
 BiocManager::install("data.table")
 BiocManager::install("dplyr")
 BiocManager::install("devtools")
 BiocManager::install("ChAMP")
 BiocManager::install("ENmix")
 BiocManager::install("CpGassoc")
 BiocManager::install("GenomicRanges")
 BiocManager::install("mCSEA")
 install_github("raivokolde/seqlm")
 install_github("WjinZhang/DMRIntTk")
 library(DMRIntTk)
 ```
 ## Running the tests
 1. Quick use of DMRIntTk.
 This is the pipeline of the idetification and integration of DMRs on 450K methylation array data:
```R
beta = load(system.file("extdata", "beta_450K.RData", package = 'DMRIntTk'))
pd = read.csv(system.file("extdata", "pd_450K.csv", package = 'DMRIntTk'))
totalDMR = identify_DMR(beta = beta, method = c("bumphunter","combp","ipDMR","mCSEA","ProbeLasso","seqlm"), pheno = pd, arraytype = "450K", group1 = "Tumor", group2 = "Normal", minProbes = 3)
bin_method = DMRInt_method(totalDMR, arraytype = "450K" )
bin_weight=DMRInt_weight(bin_method, totalDMR, pd, beta, group1 = "Tumor", group2 = "Normal")
Res=DMRInt_densitypeak(bin_weight, totalDMR, prefer = "probe", arraytype = "450K")
```
 3. Step-by-step use of DMRIntTk
### DMR sets identification
 Since DMR integration requires multiple DMR sets predicted by different methods as inputs, users should either have the self-identified multiple DMR sets, or directly use the DMR detection functions provided by DMRIntTk to identify DMR sets.
 For the latter situation, with the **methylation beta value matrix "beta"** and **the phenomenon information "pd"**, users can easily obtain the desired DMR sets with following functions(p.s.: the arraytype and minimum 
 probes can be customized, here we took 450K array and 3 probes for the example):
 ```R
 beta = load(system.file("extdata", "beta_450K.RData", package = 'DMRIntTk'))
 pd = read.csv(system.file("extdata", "pd_450K.csv", package = 'DMRIntTk'))
totalDMR = identify_DMR(beta = beta, method = c("bumphunter","combp","ipDMR","mCSEA","ProbeLasso","seqlm"), pheno = pd, arraytype = "450K", group1 = "Tumor", group2 = "Normal", minProbes = 3, regionsTypes = "promoter")
totalDMR = DMRInt_input(totalDMR, beta , group1 = "Tumor", group2 = "Normal" , arraytype = "450K")            
```

 ### DMR sets integration
 For self-identified DMR sets, users should organize them into the total DMR file containing **chromosome, start, end and methodnames**. The example total DMR file is provided.
 ```R
 beta = load(system.file("extdata","beta_450K.RData",package = 'DMRIntTk'))
 totalDMR = read.csv(system.file("extdata","totalDMR.csv",package = 'DMRIntTk'))
 ```
 #### DMRInt_input
 Prepare the input of DMRIntTk.
 This function **finds the probes included in each DMR, and calculates the methylation differences of each probe and DMR**.
 DMRInt_input function needs two files as input : 1. **Total DMR sets** obtained from different methods and 2. **Methylation level beta matrix** of all samples. 
```R
totalDMR = DMRInt_input(totalDMR, beta , group1 = "Tumor", group2 = "Normal" , arraytype = "450K")                                                        
```
 
 #### DMRInt_method
This function determines DMRs that cover each bin, and calculate the numbers of methods that cover the bin. 
It requires the **pre-split genomic bins file** and **total DMR sets** processed by the function DMRInt_input.

```R
bin_method = DMRInt_method(totalDMR, arraytype = "450K" )
```
 #### DMRInt_weight
 This function **calculates the weights of bins**. It requires **the pre-split genomic bins file**(processed from DMRIntTk_method function).
 
```R
bin_weight=DMRInt_weight(bin_method, totalDMR, pd, beta, group1 = "Tumor", group2 = "Normal")
```
 
 #### DMRInt_densitypeak
 This function **clusters all bins** based on density peak algorithm. It needs the **bins file** with calculated weights(obtained from DMRInt_weight function)
 and **total DMR sets**(obtained from DMRIntTk_input function).

```R
Res=DMRInt_densitypeak(bin_weight, totalDMR, prefer = "probe", arraytype = "450K")
```
 
 # Built With
  R is a free software environment for statistical computing and graphics.
  
 # Authors
* Xiaoqing Peng - [Central South University](https://life.csu.edu.cn/jsxx.jsp?urltype=news.NewsContentUrl&wbtreeid=1815&wbnewsid=3625)
* Wenjin Zhang - [Central South University](https://life.csu.edu.cn/)

# Version update
* **DMRIntTk** v0.1.0 -  News version 0.1.0, 2024.05.31. The first version.
