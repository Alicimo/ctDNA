library(data.table)
library(stringr)

poseidon.load.data <- function(root.dir="~/OneDrive/projects/ctDNA/",reprocess=FALSE){
  data.file <- paste0(root.dir,"data/poseidon/poseidon.Rdata")
  if(reprocess | !file.exists(data.file)){ 
    data <- poseidon.process.data(root.dir)
    save(data,file=data.file)
  } else {
    load(data.file)
  }
  data
}

poseidon.process.data <- function(root.dir="~/OneDrive/projects/ctDNA/"){

  #load sWG data
  source(paste0(root.dir,"src/load.sWG.R"))
  data <- sWG.load.data(root.dir)
  
  #Load patient meta files, subset and merge
  pat.meta <- fread(paste0(root.dir,"data/poseidon/patient.meta.csv"),check.names=T)
  
  #identify relevant samples (Emma's samples are also in set)
  x <- !is.na(match(data$meta[,Sample.name],pat.meta[,Sample.name]))
  data$cns <- data$cns[,c(rep(TRUE,4),x),with=F]
  data$segments <- data$segments[,c(rep(TRUE,4),x),with=F]
  data$meta <- data$meta[x]
  stopifnot(all.equal(data$meta[,fname],names(data$cns)[-(1:4)]))
  
  data$meta <- cbind(data$meta,pat.meta[match(data$meta[,Sample.name],pat.meta[,Sample.name]),-c("Sample.name")])
  
  print("Missing the following samples:")
  print(pat.meta[,Sample.name][!(pat.meta[,Sample.name] %in% data$meta[,Sample.name])])
  
  #easy ones - plasma + C#D#
  data$meta[type=='plasma' & grepl('C\\d*D\\d*',Sample.name),sample.time:=str_extract(Sample.name,'C\\d*D\\d*')]
  data$meta[,cycle:=as.integer(str_sub(str_extract(sample.time,'C\\d*'),2,-1))]
  data$meta[,day:=as.integer(str_sub(str_extract(sample.time,'D\\d*'),2,-1))]
  data$meta[,total.days:=(cycle-1)*28+day]
  
  #set FFPE samples to zero days
  data$meta[type=='FFPE',total.days:=0]
  
  #set all other samples to be post trail
  data$meta[,total.days:=sapply(total.days,function(x) if(is.na(x)){max(total.days,na.rm=T)+1}else{x}),patient]
  
  #either overwrite replicates with merged data or rename replicates
  if(TRUE){
    print("Overwriting duplicates with merged data")
    merged <- list()
    merged$cns <- fread(paste0(root.dir,"data/sWG/QDNA.merged/QDNA.merged.copynumbers.tsv"))
    merged$segments <- fread(paste0(root.dir,"data/sWG/QDNA.merged/QDNA.merged.segments.tsv"))
    for(sname in data$meta[,.N,Sample.name][N>1,Sample.name]){
      replicates <- data$meta[,which(Sample.name==sname)]
      
      #remove all entries bar the first
      data$meta <- data$meta[-(replicates[-1])]
      data$cns <- data$cns[,-(replicates[-1]+4),with=F]
      data$segments <- data$segments[,-(replicates[-1]+4),with=F]
      
      #overwrite remaining cns and segments with merged data
      i <- which(data$meta[replicates[1],fname] == names(merged$cns))
      data$cns[,(replicates[1]+4) := merged$cns[,i,with=F]]
      data$segments[,(replicates[1]+4) := merged$segments[,i,with=F]]
    }
    stopifnot(all.equal(data$meta[,fname],names(data$cns)[-(1:4)]))
  } else {
    print("Renaming duplicates")
    print(data$meta[,.N,Sample.name][N>1,.N])
    data$meta[,Sample.name:=if(.N>1){paste(Sample.name,1:.N,sep='-TR')}else{Sample.name},Sample.name]
  }
  
  data
}

poseidon.load.RECIST <- function(root.dir="/Users/martin06/OneDrive/projects/ctDNA/",reprocess=FALSE){
  data.file <- paste0(root.dir,"data/poseidon/poseidon.RECIST.Rdata")
  if(reprocess | !file.exists(data.file)){ 
    RECIST <- suppressWarnings(melt(fread(paste0(root.dir,"data/poseidon/patient_best_res2.csv")),id=1:2,measure=3:8,value.name="RECIST",variable.name = "cycle"))
    RECIST[,cycle:=as.integer(gsub("cyc","",cycle))]
    RECIST <- RECIST[!is.na(RECIST)]
    save(RECIST,file=data.file)
  } else {
    load(data.file)
  }
  RECIST
}