require(data.table)

sWG.load.data <- function(data.dir){
  data <- list()
  data$cns <- fread(paste0(data.dir,"QDNA.merged.copynumbers.tsv"),header=T)
  data$segments <- fread(paste0(data.dir,"QDNA.merged.segments.tsv"),header=T)

  #remove ".bam" from the fname to match the QDNA output headers
  #data$meta <- fread(paste0(data.dir,"QDNA.meta.csv"),check.names=T)
  #data$meta[,fname:=gsub('.{4}$','',fname)]
  #stopifnot(all.equal(data$meta[,fname],names(data$cns)[-(1:4)]))
 
  data
}