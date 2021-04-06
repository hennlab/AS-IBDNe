#!/usr/bin/env Rscript
# Author: Gerald Van Eeden gveeden@sun.ac.za
#Rscript msp_to_bed.R <path to msp files> </ouputdir/outputprefix> <num cores> <pop1>,<pop2>,<pop3> etc...

list.of.packages <- c("doParallel", "magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(magrittr)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

#path to MSP files
input_path = args[1]
#Output path
out_path = args[2]
#Number of cores for parallel processing
cores = args[3]
cores = strtoi(cores)
pops = args[4]
pops = unlist(strsplit(pops, ",")[[1]])

registerDoParallel(cores)

input = strsplit(input_path, split = "")
if(unlist(input)[length(input)] != "/"){
  input_path = paste(input_path, "/", sep = "")
}

files = system(paste("ls -d ",input_path,"*.msp.tsv", sep = ""), intern = T)

lapply(files, function(f){
  read.csv(f,
           header = T,
           skip = 1,
           sep = "\t")
}) ->
  msps

matrices = c()

for(m in msps){
  foreach(i=7:ncol(m))%dopar%{
    ancs = unique(m[,i])
    lapply(ancs, function(f){
      x = which(m[,i] == f)
      out = matrix(ncol = 6, nrow = nrow(m[x,]))
      out[,1:3] = unlist(m[x,1:3])
      out[,5:6] = unlist(m[x,4:5])
      out[,4] = pops[f+1]
      out
    }) %>%
      do.call("rbind", .)
  } -> msp_conv

  for(i in 1:length(msp_conv)){
    if(is.null(matrices[i][[1]])){
      matrices[[i]] = msp_conv[[i]]
    } else{
      matrices[[i]] = rbind(matrices[[i]], msp_conv[[i]])
    }
  }
}

for(i in 1:length(matrices)){
  write.table(matrices[[i]],
              paste(out_path, colnames(msps[[1]])[i+6], ".BED", sep = ""),
              quote = F,
              row.names = F,
              col.names = F)
}
