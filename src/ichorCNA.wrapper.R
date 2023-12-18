library(HMMcopy)
library(ichorCNA)
library(GenomicRanges)

normalise.ichorCNA <- function(data,indNorm){
  meds <- apply(data$cns[,indNorm+4,with=F],1,function(x){median(x,na.rm=T)})
  data$cns.normed <- cbind(data$cns[,1:4],data$cns[,.SD[,-(1:4)]] - meds)
  return(data)
}

get.env.vars.default <- function(){
  env.vars <- list(
    normal = c(0.5,0.6,0.7,0.8,0.9),
    ploidy = c(2,3),
    scStates = c(1,3),
    chrTrain = c(1:22),
    minSegmentBins = 50,
    altFracThreshold = 0.05,
    maxFracCNASubclone = 0.7,
    maxFracGenomeSubclone = 0.5,
    estimateScPrevalence = TRUE,
    lambdaScaleHyperParam = 3
  )
  return(env.vars)
}

get.env.vars.low_tumour <- function(){
  env.vars <- get.env.vars.default()
  env.vars$normal <- c(0.95, 0.99, 0.995, 0.999)
  env.vars$ploidy <- c(2)
  env.vars$estimateScPrevalence <- FALSE
  env.vars$scStates <- integer()
  return(env.vars)
}

run.ichorCNA <- function(data.normed,ind,env.vars){
  list2env(env.vars,envir = environment())
  results <- list()
  tumour_copy <- list()
  tumour_copy[["WTF"]] <- RangedData(makeGRangesFromDataFrame(data.normed$cns[,c(2:4,ind+4),with=F],keep.extra.columns = T))
  colnames(tumour_copy[[1]])[2] <- "copy"
  tumour_copy[[1]]$valid <- TRUE
  chrInd <- space(tumour_copy[[1]]) %in% chrTrain
  valid <- tumour_copy[[1]]$valid
  
  loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                                 dimnames = list(c(), c("init", "n_est", "phi_est", "Subclone_Frac", 
                                                        "Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
  logR <- as.data.frame(lapply(tumour_copy, "[[", "copy"))
  
  counter <- 1
  for(n in normal){
    for(p in ploidy){
      
      param <<- getDefaultParameters(logR[chrInd, , drop=F], ploidy = floor(p), ct.sc=scStates)
      logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
      param$lambda <- rep(logR.var, length(param$ct))
      param$lambda[param$ct %in% c(2)] <- logR.var 
      param$lambda[param$ct %in% c(1,3)] <- logR.var 
      param$lambda[param$ct >= 4] <- logR.var / 5
      param$lambda[param$ct == max(param$ct)] <- logR.var / 15
      param$lambda[param$ct.sc.status] <- logR.var / 10
      param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct))
      param$phi_0 <- p
      param$n_0 <- n
      
      hmmResults.cor <- ichorCNA::HMMsegment(tumour_copy, dataType = "copy", param = param, validInd = valid, 
                                             estimateNormal = T, estimatePloidy = T, estimateSubclone = T, verbose = F)
      
      iter <- hmmResults.cor$results$iter
      results[[counter]] <- hmmResults.cor
      loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
      subClonalBinCount <- sum(hmmResults.cor$cna[[1]]$subclone.status - 1)
      fracGenomeSub <- subClonalBinCount / nrow(hmmResults.cor$cna[[1]])
      fracAltSub <- subClonalBinCount /  sum(hmmResults.cor$cna[[1]]$copy.number != 2)
      fracAltSub <- if(is.na(fracAltSub)) {0} else {fracAltSub}
      if (subClonalBinCount == 0){
        scFrac <- NA
      }else{
        scFrac <- signif(1 - hmmResults.cor$results$sp[, iter], digits = 4)
      }
      loglik[counter, "Frac_genome_subclonal"] <- signif(fracGenomeSub, digits=2)
      loglik[counter, "Frac_CNA_subclonal"] <- signif(as.numeric(fracAltSub), digits=2)
      loglik[counter, "init"] <- paste0("n", n, "-p", p)
      loglik[counter, "n_est"] <- signif(hmmResults.cor$results$n[, iter], digits = 2)
      loglik[counter, "phi_est"] <- signif(hmmResults.cor$results$phi[, iter], digits = 4)
      loglik[counter, "Subclone_Frac"] <- scFrac
      counter <- counter + 1
    }
  }
  i <- order(as.numeric(loglik[, "loglik"]), decreasing = T)
  return(list(LLs=loglik[i,],models=results[i]))
}


check.diploid <- function(data.normed,ind,env.vars,ichorCNA.results){
  list2env(env.vars,envir = environment())

  tumour_copy <- list()
  tumour_copy[["WTF"]] <- RangedData(makeGRangesFromDataFrame(data.normed$cns[,c(2:4,ind+4),with=F],keep.extra.columns = T))
  colnames(tumour_copy[[1]])[2] <- "copy"
  tumour_copy[[1]]$valid <- TRUE
  chrInd <- space(tumour_copy[[1]]) %in% chrTrain
  valid <- tumour_copy[[1]]$valid
  
  logR <- as.data.frame(lapply(tumour_copy, "[[", "copy"))

  N <- nrow(ichorCNA.results$LLs)
  for (i in 1:N){
    hmmResults.cor <- ichorCNA.results$models[[i]]
    iter <- hmmResults.cor$results$iter
    
    ## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
    ## check if there is an altered segment that has at least a minimum # of bins
    segsS <- hmmResults.cor$results$segs[[1]]
    segsS <- segsS[segsS$chr %in% chrTrain, ]
    segAltInd <- which(segsS$event != "NEUT")
    maxBinLength = -Inf
    if (sum(segAltInd) > 0){
      maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
      maxSegRD <- RangedData(space=segsS$chr[segAltInd[maxInd]], 
                             ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
      hits <- findOverlaps(query=maxSegRD, subject=tumour_copy[[1]][valid, ])
      maxBinLength <- length(subjectHits(hits))
    }
    ## check if there are proportion of total bins altered 
    # if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
    cnaS <- hmmResults.cor$cna[[1]]
    altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
    altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
    if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
      print(paste0("Model ",i,": Failed min. req. Reseting tumour frac to 1"))
      ichorCNA.results$LLs[i,"n_est"] <- 1.0
    }
  }
  
  return(ichorCNA.results)
}

circus.plot <- function(segs, ann.chromo.frac=.25){
  
  chromo <- makeGRangesFromDataFrame(segs[[1]][,.(start=min(start),end=max(end)),chr])
  seqlengths(chromo) <- segs[[1]][,max(end),chr][,V1]
  ann.chromo <- GRanges(" ",IRanges(1, 1))
  seqlengths(ann.chromo) <- 1
  chromo.weights <- (1-ann.chromo.frac)*seqlengths(chromo)/ sum(as.numeric(seqlengths(chromo)))
  chromo.weights <- c(chromo.weights," "=ann.chromo.frac)
  chromo <- suppressWarnings(append(chromo,ann.chromo))
  
  p <- ggbio()
  
  for(i in 1:length(segs)){
    x <- makeGRangesFromDataFrame(segs[[i]],keep.extra.columns = T,seqinfo = seqinfo(chromo))
    x <- x[x$copy.number != 2]
    if(length(x)!=0){
      p <- p + circle(x,geom="rect",aes(fill=factor(copy.number),linetype=subclone.status),chr.weight=chromo.weights,trackWidth = 3,radius=8+3*i,grid=T,grid.n=1,grid.background="white",grid.line="grey70")
    } else {
      p <- p + circle(chromo,geom="segment",colour=NA,chr.weight=chromo.weights,trackWidth = 3,radius=8+3*i,grid=T,grid.n=1,grid.background="white",grid.line="grey70")
    }
    p <- p + annotate(geom="text",label=names(segs)[i],x = -1,y=9.5+3*i,hjust=1)
  }
  
  p <- p + 
    scale_fill_manual(values=c("1"="steelblue","3"="orange","4"="orangered","5"="red"),name="Copy\nNumber") + 
    scale_linetype_manual(values=c("FALSE"=1,"TRUE"=3),name="Subclone\nCNA",guide=guide_legend(override.aes = list(fill="white"))) +
    circle(chromo, geom = "scale", size = 2,chr.weight=chromo.weights,radius=8+3*(i+1)) + 
    circle(chromo, geom = "text", aes(label = seqnames), vjust = -1, size = 6,chr.weight=chromo.weights,radius=8+3*(length(segs)+3))
  p
}

run.ichor.pipeline <- function(data,
                               env.vars=get.default.env.vars.ichorCNA(),
                               ichor.dir="ichor/",
                               verbose=T){
  ###############################
  ## Run ichorCNA
  ###############################
  
  data.dir <- paste0(ichor.dir,"raw/")
  dir.create(data.dir, showWarnings = F,recursive = T)
  
  if(verbose) print("Running ichor")
  
  ichor <- list()
  ichor$raw <- lapply(1:(ncol(data)-4), function(i){
    sname <- names(data)[i+4]
    fname <- paste0(data.dir,sname,".RDS")
    if(verbose) print(paste("Working on:",basename(fname)))
    if(file.exists(fname)){
      if(verbose) print("Existing ichorCNA data found")
      return(readRDS(fname))
    }else{
      if(verbose) print("Running ichorCNA")
      x <- run.ichorCNA(data,i,env.vars)
      saveRDS(x,file = fname)
      return(x)
    }
  })
  names(ichor$raw) <-  names(data)[-(1:4)]
  
  ###############################
  ## Dediploiding
  ###############################
  
  if(verbose) print("Dedipping")
  
  data.dir <- paste0(ichor.dir,"dedip/")
  dir.create(data.dir, showWarnings = F,recursive = T)
  
  ichor$dedip <- lapply(1:length(ichor$raw), function(i){
    sname <- names(ichor$raw)[i]
    fname <- paste0(data.dir,sname,".RDS")
    if(file.exists(fname)) return(readRDS(fname))
    else{
      x <- check.diploid(data,i,env.vars,ichor$raw[[i]])
      saveRDS(x,file = fname)
      return(x)
    }
  })
  names(ichor$dedip) <- names(ichor$raw)
  
  ###############################
  ## Select the best models
  ###############################
  
  if(verbose) print("Model selection")
  
  best.models <- lapply(ichor$dedip, function(x) {
    
    if (env.vars$estimateScPrevalence){ 
      fracInd <- which(
        x$LLs[, "Frac_CNA_subclonal"] <= env.vars$maxFracCNASubclone & 
        x$LLs[, "Frac_genome_subclonal"] <= env.vars$maxFracGenomeSubclone
      )
      if (length(fracInd) > 0){
        return( fracInd[order(x$LLs[fracInd, "loglik"], decreasing=TRUE)][1] )
      }else{
        return(1)
      }
    }else{
      return(1)
    }
  })
  
  ichor$summary <- rbindlist(lapply(1:length(ichor$dedip),function(i) ichor$dedip[[i]]$LLs[best.models[[i]],]))
  ichor$summary$Sample.name <- names(ichor$dedip)
  
  ichor$bestModels <- lapply(1:length(ichor$dedip),function(i) ichor$dedip[[i]]$models[[best.models[[i]]]])
  names(ichor$bestModels) <- names(ichor$dedip)
  
  ichor$raw <- NULL
  ichor$dedip <- NULL
  return(ichor)
}