
# DMRIntTk
 DMRIntTk is a toolkit for integrating DMR sets predicted by different methods on a same methylation array dataset based on density peak clustering algorithm.
 It contains five main functions including DMRIntTk_input, DMRIntTk_matrix, DMRIntTk_method, DMRIntTk_weight and DMRIntTk_densitypeak.
The following tutorial sections will teach users to use the DMRIntTk package efficiently and accurately.
 
 ##Getting started 
 
 These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
 ###Prerequisites
 You need to install R (version >= 2.10).

 ##Installing
 You can input the following commands to install the DMRIntTk R package.
 
 ```R
 install.packages("dplyr")
 install.packages("devtools")
 library(dplyr)
 library(devtools)
 install_github("WjinZhang/DMRIntTk")
 ```
 ##Running the tests
 ###DMRInt_input
 Prepare the input of DMRIntTk.
 This function finds the probes included in each DMR, and calculates the methylation differences of each probe and DMR.
 DMRIntTk_input function needs two files as input : 1. Total DMR sets  obtained from different methods and 2. Methylation level beta matrix of all samples.
 
 
```R
totalDMR = read.csv(system.file("extdata","totalDMR.csv",package = 'DMRIntTk'))
beta = readRDS(system.file("extdata","beta.RDS",package = 'DMRIntTk'))
case = 1:5
control = 6:10
totalDMR=DMRInt_input(totalDMR, beta, case,control, arraytype = "450K")
```
 
 ###DMRInt_matrix
 This function calculates the DMRscore and weights of all methods under different methylation difference thresholds and
 constructs the reliability weight matrix for each method. It requires the DMR sets file processed by the function DMRInt_input.

 ```R
 weight_m = DMRInt_matrix(totalDMR)
 ```
 
 ###DMRInt_method
This function determines DMRs that cover each bin, and calculate the numbers of methods that cover the bin. 
It requires the  pre-split genomic bins file and the DMR sets file processed by the function DMRInt_input.

```R
interval_method = DMRInt_method(totalDMR, arraytype = "450K" )
```
 ###DMRInt_weight
 This function calculates the weights of bins. It requires the pre-split genomic bins file(processed from DMRIntTk_method function) 
 and the weight matrix obtained from DMRInt_matrix function.
 
 
```R
interval_weight=DMRInt_weight(interval_method,weight_m,beta,case,control)
```
 
 
 ###DMRInt_densitypeak
 This function clusters all bins based on density peak algorithm. It needs the bins file with calculated weights(obtained from DMRInt_weight function)
 and the DMR sets file(obtained from DMRIntTk_input function).

```R
Res=DMRInt_densitypeak(interval_weight, totalDMR, prefer = "probe", arraytype = "450K")
```
 
 ##Built With
  R is a free software environment for statistical computing and graphics.
  
 ##Authors
Wenjin Zhang - Thesis code writing work Central South University
Edit By [MaHua](http://mahua.jser.me)
