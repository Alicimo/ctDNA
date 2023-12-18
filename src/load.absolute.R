load.absolute.data <- function(project.dir,meta){
  
  absolute <- list()
  absolute$meta <- fread(paste0(project.dir,'data/ABSOLUTE/reviewed/summary.martin06.ABSOLUTE.table.txt'),check.names = T)
  absolute$meta <- cbind(absolute$meta,meta[match(absolute$meta[,sample],gsub('/','-',meta[,Sample.name]))])
  absolute$meta <- absolute$meta[order(patient,sample)]
  
  absolute$samples <- list()
  for(sname in absolute$meta[,sample]){
    absolute$samples[[sname]] <- fread(paste0(project.dir,'data/ABSOLUTE/reviewed/SEG_MAF/',sname,'.segtab.txt'))
  }
  
  absolute$patients <- list()
  for(pat in unique(absolute$meta[,patient])){
    x <- rbindlist(absolute$samples[absolute$meta[,which(patient==pat)]])
    x$ploidy <- absolute$meta[match(x[,sample],absolute$meta[,sample]),ploidy]
    x[,c("gain","loss"):=.(if(ploidy>2.7){modal_cn>8}else{modal_cn>4},if(ploidy>2.7){modal_cn<(ploidy-2.7)}else{modal_cn==0}),.(sample,Chromosome,Start.bp)]
    absolute$patients[[pat]] <- x
  }
  
  absolute
}