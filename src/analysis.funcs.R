construct.survival.data <- function(detect){
  
  x <- detect$RECIST[RECIST=="PD",head(.SD,1),Patient.ID]
  x <- rbind(x,detect$RECIST[!(Patient.ID %in% x$Patient.ID),tail(.SD,1),Patient.ID])
  x <- x[,.(Patient.ID,event.date=date,RECIST)]
  x <- x[!is.na(RECIST)]
  
  y <- get.measures.meta(detect)
  y <- y[Sample.type=="ctDNA"][order(Patient.ID,Sample.num)]
  y[,Date.collected.end:=c(tail(Date.collected,-1),NA),Patient.ID]
  
  x <- merge(x,y,by="Patient.ID")
  
  #Fix patients where the first event (PD) is before any samples.
  # In this scenario, remove that event from being considered and look for a future event
  p.IDs <- x[,floor(difftime(event.date,min(Date.collected),units="days")),Patient.ID][V1 < 1,unique(Patient.ID)]
  old.p.IDs <- x[,unique(Patient.ID)]
  while(!all(suppressWarnings(old.p.IDs==p.IDs))){
    y <- detect$RECIST[Patient.ID %in% p.IDs][RECIST=="PD",.SD[2,],Patient.ID][!is.na(date)]
    y <- rbind(y,detect$RECIST[Patient.ID %in% setdiff(p.IDs,y$Patient.ID),tail(.SD,1),Patient.ID])
    
    for(p.ID in p.IDs){
      x[Patient.ID==p.ID,event.date:=y[Patient.ID==p.ID,date]]
      x[Patient.ID==p.ID,RECIST:=y[Patient.ID==p.ID,RECIST]]
    }
    
    old.p.IDs <- p.IDs
    p.IDs <- x[,floor(difftime(event.date,min(Date.collected),units="days")),Patient.ID][V1 < 1,unique(Patient.ID)]
  }
  
  x <- x[,.(
    start=floor(difftime(Date.collected,min(Date.collected),units="days")),
    end=floor(difftime(Date.collected.end,min(Date.collected),units="days")),
    event.date=floor(difftime(event.date,min(Date.collected),units="days")),
    ichorCNA.best,ichorCNA.default,ichorCNA.low,ichorCNA.suggested,z.OR,t.Mad,absolute,RECIST),
    Patient.ID]
  
  # Add genotype
  x <- merge(x,detect$genotype,by="Patient.ID")
  x[,ER.Her2:=interaction(ER.status,Her2.status)]
  
  # Remove patients that progress before plasma samples
  x <- x[event.date>0]
  
  # Remove samples taken after PD
  x <- x[event.date>start]
  
  # Set final sample to "finish" on final scan (if na)
  x[is.na(end),end:=event.date]
  
  # Set events
  x[,event:=0]
  x[(event.date>start & event.date <= end) & RECIST == "PD", event:=1]
  
  #If final ctDNA sample extends past a PD call, set end to event.date
  x[event.date<end & event==1,end:=event.date]
  
  return(x)
}

#################################################################################

time.dependent.conc <- function(model,dt,id.name){
  strata <- dt[[id.name]]
  x <- survfit(model,newdata=dt,id=strata)
  y <- dcast(data.table(
    id=as.numeric(rep(names(x$strata),x$strata)),
    time=x$time,
    cumhaz=x$cumhaz
  ),time~id)
  z <- dt[,.(sum(event),tail(event.date,1)),id.name]
  
  n <- length(z$patient.ID)
  conc <- as.list(table(unlist(lapply(1:n, function(i){
    sapply(i:n,function(j){
      if(i==j) return(NA)
      k <- c(i,j)
      if(z[k,sum(V1)]==0) return(NA)
      if(z[k,sum(V1)]==1 &  z[k,order(V1) == order(V2)][1]) return(NA)
      if(z[k,sum(V1)]==1 &  (z[k[1],V2]==z[k[2],V2])) return(NA) #This does not have to be NA
      
      x <- cbind(
        z[k],
        V3=as.numeric(y[which(apply(is.na(y[,k+1,with=F]),1,any))[1]-1,k+1,with=F])
      )
      
      if(x[1,V3]==x[2,V3]) return("Tied")
      return(ifelse(x[,order(V2)==order(V3)][1],"Disagree","Agree"))
    })
  }))))
  return((conc$Agree + .5*conc$Tied) / sum(unlist(conc)))
}

#################################################################################

get.treat.bound <- function(detect,post.padding=0){
  x <- detect$treatment
  
  pre <- rbind(
    x[,.(treatment="unknown",start.date=as.Date("0001-01-01"),end.date=start.date[1]),Patient.ID],
    x[,head(.SD,-1),Patient.ID]
  )[order(Patient.ID,start.date)]
  post <- x
  
  names(pre) <- paste0("pre.",names(pre))
  names(post) <- paste0("post.",names(post))
  
  x <- cbind(pre,post)
  x <- x[,.(patient.id=pre.Patient.ID,pre.treatment,post.treatment,
            pre.start.date=pre.start.date,
            change.date=pre.end.date,
            post.end.date=post.end.date)]
  
  x <- rbindlist(lapply(1:nrow(x),function(i){
    tb <- x[i]
    
    pre.score <- merge(
      detect$sWG$meta[Sample.type=="ctDNA"][Patient.ID==tb$patient.id & (Date.collected > tb$pre.start.date & Date.collected < tb$change.date)][order(Date.collected)][.N],
      detect$ichorCNA.default$summary
    )[,.(pre.score.date=Date.collected-tb$change.date,pre.score=1-n_est)][1]
    
    post.score <- merge(
      detect$sWG$meta[Sample.type=="ctDNA"][Patient.ID==tb$patient.id & (Date.collected < tb$post.end.date & Date.collected > (tb$change.date+post.padding))][order(Date.collected)][1],
      detect$ichorCNA.default$summary
    )[,.(post.score.date=Date.collected-tb$change.date,post.score=1-n_est)][1]
    
    pre.CT <- detect$RECIST[Patient.ID==tb$patient.id & (date > tb$pre.start.date & date < tb$change.date)][order(date)][.N][,.(pre.CT.date=date-tb$change.date,pre.CT=RECIST)][1]
    
    post.CT <- detect$RECIST[Patient.ID==tb$patient.id & (date < tb$post.end.date & date > (tb$change.date+post.padding))][order(date)][1][,.(post.CT.date=date-tb$change.date,post.CT=RECIST)][1]
    
    y <- list(tb,pre.score,post.score,pre.CT,post.CT)
    return(do.call("cbind",y))
  }))
  x[,dif.score:=post.score-pre.score]
  return(x)
}

get.treat.bound.preprocess <- function(detect, post.padding=0, merge.CT=FALSE, m.terms=c(), time.limit=-1){
  treat.bound <- get.treat.bound(detect,post.padding)

  if(merge.CT){
    treat.bound[,pre.CT:=factor(ifelse(pre.CT=="PD","PD","!PD"))]
    treat.bound[,post.CT:=factor(ifelse(post.CT=="PD","PD","!PD"))]
  }

  if(length(m.terms)>0){
    treat.bound <- treat.bound[,names(treat.bound) %in% m.terms,with=F]
    treat.bound <- treat.bound[rowSums(is.na(treat.bound))==0]
  }
  
  treat.bound <- treat.bound[rowSums(treat.bound[,grep("(CT|score).date",names(treat.bound)),with=F][,lapply(.SD,abs)] > time.limit) == 0]
  
  return(treat.bound)
}

#################################################################################

get.dist.mat <- function(segs.expanded,project.dir,dist.func,...){
  fname <- paste0(project.dir,"data/detect.dist_mat.RDS")
  if(file.exists(fname)){
    dist.matrix <- readRDS(fname)
  } else {
    x <- segs.expanded[,-(1:2)]
    N <- ncol(x)
    dist.matrix <- matrix(0,N,N)
    colnames(dist.matrix) <- colnames(x)
    rownames(dist.matrix) <- colnames(x)
    for(i in 1:N){
      for(j in i:N){
        d <- dist.func(x[,i,with=F],x[,j,with=F],...)
        dist.matrix[i,j] <- d
        dist.matrix[j,i] <- d
      }
    }
    saveRDS(dist.matrix,fname)
  }
  return(dist.matrix)
}

cor.dist <- function(x,y,...){
  if(var(x,na.rm=T)==0 | var(y,na.rm=T)==0) return(NA)
  else return(cor(x,y,...))
}

custom.dist <- function(x,y,...){
  x <- sum(x != y,na.rm = T)
  x <- x - 0.5*sum( (x>2 & y>2) & (x != y) ,na.rm = T)
  x <- x / sum( !is.na(x) & !is.na(y) )
}

#################################################################################

expand.segmentation <- function(segs){
  segs.expanded <- lapply(segs, function(x){
    x[,.(seq(start-1,end-1e6,1E6),copy.number),.(chr,start,end)][,.(chr,start=V1,copy.number)]
  })
  if(length(segs.expanded) == 1) segs.expanded <- segs.expanded[[1]]
  else {
    i <- which.max(sapply(segs.expanded,nrow))
    segs.expanded <- cbind(
      segs.expanded[[i]][,-3],
      sapply(segs.expanded, function(x){
        merge(segs.expanded[[i]],x,by=c("chr","start"),all.x=TRUE,sort=F)[,copy.number.y]
      })
    )
  }
  return(segs.expanded)
}

aggregate.segmentation <- function(segs.expanded){
  x <- segs.expanded[,-(1:2)]
  N <- ncol(x)
  x <- cbind(segs.expanded[,1:2],V1=rowSums(x<2,na.rm=T),V2=rowSums(x>2,na.rm=T),V3=rowSums(is.na(x)))
  
  x[V3>(N/10),V1:=NA]
  x[,V1:=V1/(N-V3)]
  x[,V1:=zoo::na.locf(V1)]
  
  x[V3>(N/10),V2:=NA]
  x[,V2:=V2/(N-V3)]
  x[,V2:=zoo::na.locf(V2)]
  
  x
}

reduce.segmentation.wrapper <- function(segs.aggregated){
  l <- segs.aggregated[,.(chr,start,V1)]
  g <- segs.aggregated[,.(chr,start,V1=V2)]
  x <- list(
    gain =  reduce.segmentation(g),
    loss =  reduce.segmentation(l)
  )
  x
}

reduce.segmentation <- function(x){
  x$grp <- x[,c(1,which(diff(V1)!=0)+1,.N+1),chr][,rep(1:length(diff(V1)), diff(V1)),chr][,V1]
  x <- x[,.(start=min(start)+1,end=max(start)+1E6,V1=unique(V1),N.bins=.N),.(chr,grp)][,-c('grp')]
  x$chr <- factor(x$chr,levels = c(1:22,"X"))
  x
}

process.segs <- function(segs.raw){
  segs <- list()
  segs$raw <- segs.raw
  segs$expanded <- expand.segmentation(segs$raw)
  segs$aggregated <- aggregate.segmentation(segs$expanded)
  segs$agg.reduced <- reduce.segmentation.wrapper(segs$aggregated)
  return(segs)
}

detect.extract.segs <- function(detect){
  segs <- lapply(seq_along(detect$ichorCNA.default$bestModels), function(i){
    x <- data.table(detect$ichorCNA.default$bestModels[[i]]$results$segs[[1]])
    if(detect$ichorCNA.default$summary[i,n_est]==1) x[,copy.number:=2]
    return(x)
  })
  names(segs) <- names(detect$ichorCNA.default$bestModels)
  return(segs)
}

get.detect.best_ctDNA_index <- function(detect){
  x <- cbind(detect$sWG$meta,detect$ichorCNA.default$summary)
  i <- x[,.I[Sample.type=="ctDNA" & n_est!=1]]
  i <- i[ x[i][,.I[which.min(n_est)],Patient.ID][,V1] ]
  #segs.ctDNA <- process.segs(detect$ichorCNA.default$bestModels$raw[i])
  return(i)
}

#################################################################################

generate.her2.agg <- function(her2.samples,gain.thresh=0.4,loss.thresh=-0.5){
  
  her2.samples <- by(her2.samples,her2.samples[,sample],function(x){
    x <- x[,.(chr=gsub("chr","",chr),start=round(start/1E6)*1E6 + 1,end=round(end/1E6)*1E6,copy.number=log2((nMajor + nMinor)/ploidy))]

    x <- x[start<end]
    if(x[,diff(end)>=0,chr][V1==F,.N]) print("HERE")
    if(x[,tail(start,-1)-head(end,-1),chr][V1<0,.N]) print("Also HERE")
    x[,copy.number:=ifelse(copy.number>gain.thresh,3,ifelse(copy.number<loss.thresh,1,2))]
    x[chr==23,chr:="X"]
    data.table(x)
  })
  
  her2.samples <- process.segs(her2.samples)
  return(her2.samples)
}

trim.segmentation <- function(segs,chr.max,chr.min){
  for(i in 1:nrow(chr.max)){
    segs$agg.reduced$gain <- segs$agg.reduced$gain[(chr == chr.max[i,chr] & start < chr.max[i,V1]) | chr != chr.max[i,chr]]
    segs$agg.reduced$gain[chr == chr.max[i,chr] & end > chr.max[i,V1],end:=chr.max[i,V1]]
    segs$agg.reduced$loss <- segs$agg.reduced$loss[(chr == chr.max[i,chr] & start < chr.max[i,V1]) | chr != chr.max[i,chr]]
    segs$agg.reduced$loss[chr == chr.max[i,chr] & end > chr.max[i,V1],end:=chr.max[i,V1]]
  }
  
  for(i in 1:nrow(chr.min)){
    segs$agg.reduced$gain <- segs$agg.reduced$gain[(chr == chr.min[i,chr] & end > chr.min[i,V1]) | chr != chr.min[i,chr]]
    segs$agg.reduced$gain[chr == chr.min[i,chr] & start < chr.min[i,V1],start:=chr.min[i,V1]]
    segs$agg.reduced$loss <- segs$agg.reduced$loss[(chr == chr.min[i,chr] & end > chr.min[i,V1]) | chr != chr.min[i,chr]]
    segs$agg.reduced$loss[chr == chr.min[i,chr] & start < chr.min[i,V1],start:=chr.min[i,V1]]
  }
  return(segs)
}

calc.gain.frac.diff <- function(gain.thresh,her2.samples,chr.max,chr.min,total.gain){
  x <- generate.her2.agg(her2.samples,gain.thresh=gain.thresh)
  x <- trim.segmentation(x,chr.max,chr.min)
  return(abs(x$agg.reduced$gain[,sum((end-start)*V1)] - total.gain))
}

calc.loss.frac.diff <- function(loss.thresh,her2.samples,chr.max,chr.min,total.loss){
  x <- generate.her2.agg(her2.samples,loss.thresh = loss.thresh)
  x <- trim.segmentation(x,chr.max,chr.min)
  return(abs(x$agg.reduced$loss[,sum((end-start)*V1)] - total.loss))
}

get.metabric <- function(root.dir,project.dir,detect){
  fname <- paste0(project.dir,"data/metabric.RDS")
  if(file.exists(fname)){
    metabric <- readRDS(fname)
  } else {
    
    metabric <- list(
      meta = fread(paste0(root.dir,"metabric/data/patientData.txt")),
      cna = fread(paste0(root.dir,"metabric/data/ascatSegments.txt"))
    )
    
    if(1){
      metabric$segs <- generate.her2.agg(metabric$cna)
      metabric$opt.thr <- list(gain=0.4,loss=-0.5)
    } else {
      ctDNA.models <- process.segs(detect$ichorCNA.default$bestModels$raw[get.detect.best_ctDNA_index(detect)])
      
      chr.max <- ctDNA.models$agg.reduced$gain[,max(end),chr]
      chr.min <- ctDNA.models$agg.reduced$gain[,min(start),chr]
  
      total.gain <- ctDNA.models$agg.reduced$gain[,sum((end-start)*V1)]
      total.loss <- ctDNA.models$agg.reduced$loss[,sum((end-start)*V1)]
      
      opt.thr.gain <- optimise(calc.gain.frac.diff,metabric$cna,chr.max,chr.min,total.gain,interval = c(0,1),tol=1E-3)
      opt.thr.loss <- optimise(calc.loss.frac.diff,metabric$cna,chr.max,chr.min,total.loss,interval = c(-1,0),tol=1E-3)
    
      metabric$segs <- generate.her2.agg(metabric$cna,
                                         opt.thr.gain[[1]],
                                         opt.thr.gain[[2]])
      metabric$opt.thr = list(
        gain = opt.thr.gain[[1]],
        loss = opt.thr.loss[[1]]
      )
    }
    saveRDS(metabric,fname)
  }
  return(metabric)
}

#################################################################################

get.gl.counts <- function(segs.expanded){
  x <- segs.expanded[,-(1:2)]
  N <- ncol(x)
  x <- cbind(segs.expanded[,1:2],V1=rowSums(x<2,na.rm=T),V2=rowSums(x>2,na.rm=T),V3=rowSums(is.na(x)))
  
  x[V3>(N/10),V1:=NA]
  x[V3>(N/10),V2:=NA]
  x[is.na(V2),V3:=NA]
  
  x[,V1:=zoo::na.locf(V1)]
  x[,V2:=zoo::na.locf(V2)]
  x[,V3:=zoo::na.locf(V3)]
  
  x[,V3:=N-V3]
  
  return(x)
}
