# sWG reformat data - Script for Meiling
# reformats the QDNA data structure so that only the start and end point of segments are given
# author: Alistair Martin
# date: 1st Sep 2017

library(data.table)
load("poseidon.Rdata")

x <- lapply(1:nrow(data$meta),function(i){
  
  #extract a single sample
  x <- data$segments[,c(2:4,i+4),with=F]
  names(x)[4] <- c("segment.value")
  
  #find points where the segment value changes within each chromosome
  x$grp <- x[,c(1,which(diff(segment.value)!=0)+1,.N+1),chromosome][,rep(1:length(diff(V1)), diff(V1)),chromosome][,V1]
  
  #reduce data structure to give only the start and end point of segments
  x <- x[,.(start=min(start),end=max(end),segment.value=unique(segment.value),N.bins=.N),.(chromosome,grp)][,-c('grp')]
  x  
})