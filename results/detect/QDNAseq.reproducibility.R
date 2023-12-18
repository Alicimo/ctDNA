log2adhoc <- function(x, offset=.Machine$double.xmin, inv=FALSE) {
  if (!inv) {
    x[x < 0] <- 0
    x <- x + offset
    log2(x)
  } else {
    x <- 2^x
    x - offset
  }
}

postsegnormalize <- function(seg,inter=c(-0.1,0.1)){

  values <- c()
  for (i in 1:ncol(seg)) {
    values <- c(values, median(seg[,i]));
  }    
  matrixValues    <- matrix(rep(values, nrow(seg)), ncol=ncol(seg), byrow=TRUE);
  seg <- seg - matrixValues #postseg works best when data are median normalized 
  countlevall <- apply(seg,2,function(x) {as.data.frame(table(x))})
  
  intcount <- function(int,sv){
    sv1 <- as.numeric(as.vector(sv[,1]))
    wh <- which(sv1<=int[2] & sv1>=int[1])
    return(sum(sv[wh,2]))
  }
  
  postsegnorm <- function(segvec,int=inter,intnr=3){
    intlength <- (int[2]-int[1])/2
    gri <- intlength/intnr
    intst <- int[1]+(0:intnr)*gri
    intend <- intst+intlength
    ints <- cbind(intst,intend)
    intct <- apply(ints,1,intcount,sv=segvec)
    whmax <- which.max(intct)
    return(ints[whmax,]) 
  }
  
  postsegnorm_rec <- function(segvec,int,intnr=3){
    newint <- postsegnorm(segvec,int,intnr)
    newint <- postsegnorm(segvec,newint,intnr)
    newint <- postsegnorm(segvec,newint,intnr)
    newint <- postsegnorm(segvec,newint,intnr)
    newint <- postsegnorm(segvec,newint,intnr)
    return(newint[1]+(newint[2]-newint[1])/2)
  }
  listres <- lapply(countlevall,postsegnorm_rec,int=inter)
  vecres <- c();for(i in 1:length(listres)){vecres <- c(vecres,listres[[i]])}
  
  seg.new <- t(t(seg)-vecres)
  #cn.new <- t(t(cn-matrixValues)-vecres)
  return(seg)
}

segment.wrapper <- function(x,seed=NULL,...){
  if(!is.null(seed)) set.seed(seed)
  y <- CNA(
    genomdat=x[,5,drop=F],
    chrom=factor(x$chromosome,levels=unique(x$chromosome),ordered=T),
    maploc=x$start,
    data.type = "logratio",
    sampleid=colnames(x)[5],
    presorted = T
  )
  #y <- segment(y, alpha=1e-10, undo.splits="sdundo",undo.SD=1.0, verbose=0,...)
  y <- segment(y,...)
  seg <- as.matrix(rep(y$output$seg.mean,y$output$num.mark))
  seg <- postsegnormalize(seg)
  y <- cbind(x[,1:4],"segments"=as.vector(seg))
  y$chromosome <- factor(y$chromosome,unique(y$chromosome),ordered=T)
  return(y)
}

plot.segs <- function(log2r,segs){
  if(names(log2r)[5] != "cns") names(log2r)[5] <- "cns"
  log2r$chromosome <- factor(log2r$chromosome,unique(log2r$chromosome),ordered=T)
  p <- ggplot(x) + geom_point(data=log2r,aes(x=start,y=cns),alpha=0.5,size=0.5)
  cols <- brewer.pal(length(segs),"Set1")
  for(i in 1:length(segs)) p <- p + geom_step(data=segs[[i]],aes(x=start,y=segments),size=0.5,col=cols[i])
  p <- p + facet_wrap(~chromosome,scales="free_x") + theme_bw()
  return(p)
}

root.dir <- "~/OneDrive/projects/"
project.dir <- paste0(root.dir,"ctDNA/")

source(paste0(project.dir,"src/load.detect.R"))
source(paste0(project.dir,"src/ichorCNA.wrapper.R"))

pkgs <- c("data.table","ggplot2","ggthemes","QDNAseq","DNAcopy","Rcolorbrewer")
pkgs <- lapply(pkgs,function(pkg) suppressWarnings(suppressMessages(require(pkg, ch = TRUE))))

## FAKE DATA

run.segment.test <- function(bins,noise,cellularity,n,...){
  set.seed(42)
  genomdat <- rnorm(bins, sd=noise) + cellularity * rep(c(-1,0,1),c(floor(bins/3),floor(bins/3),floor(bins/3)+(bins%%3)))
  chrom <- rep(1,c(bins))
  maploc <- c(1:bins)
  set.seed(Sys.time())
  x <- sum(sapply(1:n,function(i) nrow(segment(CNA(genomdat, chrom, maploc),...)$output))!=0)
  return(x)
}

lapply(50*(1:4),function(n) sapply(0.1*(0:10),function(stn) run.segment.test(n,stn,1,100,verbose=0,alpha=1e-10, undo.splits="sdundo",undo.SD=1.0)))

## DETECT data

detect <- detect.load(project.dir)
detect$sWG.normed <- normalise.ichorCNA(detect$sWG,which(detect$sWG.meta[,Sample.type=="BC"]))

detect$sWG.normed$cns[,lapply(.SD[,-(1:4)],function(x) median(abs(unlist(by(x,chromosome,diff)))))]


results <- lapply(10000*(1:1),function(n){
  x <- lapply(1:nrow(detect$sWG.meta),function(i){ 
    lapply(1:2,function(j){
      segment.wrapper(detect$sWG$cns[,c(1:4,i+4),with=F],alpha=1e-4, undo.splits="sdundo",undo.SD=1.0, verbose=0,nperm = n)
    })
  })
  y <- lapply(x, function(z) lapply(z,function(y) y[,sum(diff(segments)!=0),chromosome]))
  z <- do.call(cbind,lapply(y, function(z) z[[1]][,2] - z[[2]][,2] ))
  return(z)
})
lapply(results, function(z) table(colSums(abs(z))))
sort(rowSums(z==0)/ncol(z))
table(colSums(abs(z)))

xtable(x)
x <- cbind(detect$sWG.meta,x)
ggplot(x) + aes(x=x,fill=Sample.type) + geom_density(alpha=0.5,adjust=2)
ggplot(x) + aes(x=factor(x),y=Mapped.reads) + geom_boxplot() + facet_grid(Sample.type~.) + scale_y_log10()

x$noise <- detect$sWG$cns[,sapply(.SD[,-(1:4)],function(y){median(abs(diff(y)),na.rm=T)})]
ggplot(x) + aes(x=factor(x),y=noise) + geom_boxplot() + facet_grid(Sample.type~.)

i <- 343 #lowest n_est
x <- segment.wrapper(detect$sWG$cns[,c(1:4,i+4),with=F])
y <- segment.wrapper(detect$sWG.normed$cns[,c(1:4,i+4),with=F])
x[,which(diff(segments)!=0),chromosome]
y[,which(diff(segments)!=0),chromosome]

i <- detect$sWG.meta[,which.max(Mapped.reads)]
x <- segment.wrapper(detect$sWG$cns[,c(1:4,i+4),with=F])
y <- segment.wrapper(detect$sWG.normed$cns[,c(1:4,i+4),with=F])
x[,which(diff(segments)!=0),chromosome]
y[,which(diff(segments)!=0),chromosome]


p <- lapply(10^(3:5), function(n) plot.segs(detect$sWG$cns[,c(1:4,i+4),with=F],lapply(1:9,function(j) segment.wrapper(detect$sWG$cns[,c(1:4,i+4),with=F],nperm=n))))
p


