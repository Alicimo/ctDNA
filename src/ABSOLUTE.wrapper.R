require(ABSOLUTE)

run.abs <- function(x,sample.name='unknown.sample',absolute.dir="ABSOLUTE/"){
  genome <- "hg18"
  platform <- "Illumina_WES"
  primary.disease <- 'CANCER'
  sigma.p <- 0
  max.sigma.h <- 0.02
  min.ploidy <- 0.95
  max.ploidy <- 10
  max.as.seg.count <- 1E10
  max.non.clonal <- 100
  max.neg.genome <- 0
  copy_num_type <- "total"
  results.dir <- paste0(absolute.dir,'indv/')
  tmp.fname <- tempfile()
  
  dir.create(results.dir,showWarnings=F,recursive = T)
  write.table(x,file=tmp.fname,quote=F,row.names=F,sep='\t')
  RunAbsolute(tmp.fname, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, 
              platform, sample.name, results.dir, max.as.seg.count, max.non.clonal, 
              max.neg.genome, copy_num_type, verbose=TRUE)
}

run.absolute.pipeline <- function(data,absolute.dir="ABSOLUTE/"){
  
  for(i in 1:(ncol(data)-4)){
    j <- i + 4
    sname <- colnames(data)[j]
    if( !file.exists(paste0(absolute.dir,'indv/',sname,'.ABSOLUTE.RData')) ){
      print(sname)
      x <- data[,c(2:4,j),with=F]
      names(x) <- c("Chromosome","Start","End","Segment_Mean")
      x$grp <- x[,c(1,which(diff(Segment_Mean)!=0)+1,.N+1),Chromosome][,rep(1:length(diff(V1)), diff(V1)),Chromosome][,V1]
      x <- x[,.(Start=min(Start),End=max(End),Segment_Mean=unique(Segment_Mean),Num_Probes=.N),.(Chromosome,grp)][,-c('grp')]
      x[Chromosome=="X",Chromosome:="23"]
      capture.output(invisible(run.abs(x,sname)),file = NULL)
    }
  }
  
  abs.files <- list.files(paste0(absolute.dir,"indv/"),"*.ABSOLUTE.RData",full.names=T)
  time.stamps <- lapply(abs.files, function(x) as.Date(file.info(x)$mtime))
  time.stamp.newest <- time.stamps[[which.max(time.stamps)]]
  
  review.fname <- paste0(absolute.dir,"reviews/reviews.PP-calls_tab.txt")
  time.stamp.review <- as.Date(file.info(review.fname)$mtime)
  
  if(time.stamp.newest > time.stamp.review) {
    cat("New data objects detected\nCreating new review object")
    CreateReviewObject('reviews',abs.files,paste0(absolute.dir,"reviews/"),'total',verbose=TRUE)
    ExtractReviewedResults(review.fname,'martin06',paste0(absolute.dir,"reviews/reviews.PP-modes.data.RData"),'ABSOLUTE','summary','total',verbose=TRUE)
  }
  
  x <- fread(paste0(absolute.dir,"reviews/reviews.PP-calls_tab.txt"),check.names = T)
  names(x)[1] <- "Sample.name"
  x[call.status=="low purity",purity:=0]
  return(x)
}
