source(paste0(project.dir,"src/load.detect.R"))             
source(paste0(project.dir,"src/ichorCNA.wrapper.R"))
source(paste0(project.dir,"src/ABSOLUTE.wrapper.R"))

detect.load.measures <- function(project.dir){
  fname <- paste0(project.dir,"data/detect.measures.RDS")
  if(file.exists(fname)){
    detect <- readRDS(fname)
  } else {
    detect <- detect.load(project.dir)
    detect$sWG <- normalise.ichorCNA(detect$sWG,which(detect$sWG$meta[,Sample.type=="BC"]))
    
    ### Run IchorCNA (default)
    ichor.dir <- paste0(project.dir,"data/sequencing/ichor/ichor.default/")
    env.vars <- get.env.vars.default()
    detect$ichorCNA.default <- run.ichor.pipeline(detect$sWG$cns.normed, env.vars, ichor.dir, verbose=T)
    detect$ichorCNA.default$summary <- merge(detect$sWG$meta,detect$ichorCNA.default$summary,sort=F)
    
    ### Run IchorCNA (low tumour content)
    ichor.dir <- paste0(project.dir,"data/sequencing/ichor/ichor.low/")
    env.vars <- get.env.vars.low_tumour()
    detect$ichorCNA.low <- run.ichor.pipeline(detect$sWG$cns.normed, env.vars, ichor.dir, verbose = T)
    detect$ichorCNA.low$summary <- merge(detect$sWG$meta,detect$ichorCNA.low$summary,sort=F)
    
    ### Take best models (by LL) of both ichorCNA modes
    x <- lapply(1:nrow(detect$ichorCNA.default$summary),function(i){
      if(detect$ichorCNA.default$summary[i,loglik] > detect$ichorCNA.low$summary[i,loglik]){
        return(detect$ichorCNA.default$summary[i])
      } else {
        return(detect$ichorCNA.low$summary[i])
      }
    })
    detect$ichorCNA.optimal <- list("summary" = rbindlist(x))
  
    ### Replace ichorCNA default with ichorCNA low if former predicts <5%
    x <- lapply(1:nrow(detect$ichorCNA.default$summary),function(i){
      if(detect$ichorCNA.default$summary[i,n_est] <= 0.95){
        return(detect$ichorCNA.default$summary[i])
      } else {
        return(detect$ichorCNA.low$summary[i])
      }
    })
    detect$ichorCNA.suggested <- list("summary" = rbindlist(x))
    
    ### Run ABSOLUTE
    absolute.dir <- paste0(project.dir,'data/sequencing/absolute/')
    dir.create(absolute.dir,showWarnings=F)
    x <- run.absolute.pipeline(detect$sWG$segments,absolute.dir)
    detect$absolute <- merge(detect$sWG$meta,x,all.x=T,sort=F)
    
    ### Calculate Z.OR Score
    #https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0041162
    #http://clinchem.aaccjnls.org/content/61/6/838#ref-19
    #https://onlinelibrary.wiley.com/doi/pdf/10.1002/ijc.31397
    
    x <-  cbind(detect$sWG$cns.normed[,1:4,with=F],apply(detect$sWG$cns.normed[,-(1:4)],2,function(y){y[y < -2] <- NA; y}))
    i <- which(detect$sWG$meta[,Sample.type=="BC"])
    s.OR <- apply(x[,i+4,with=F],1,function(y){sd(y,na.rm=T)})
    u.OR  <- rowSums(x[,i+4,with=F],na.rm=T) / length(i)
    detect$z.OR <- cbind(
      detect$sWG$meta,
      z.OR=x[,sapply(.SD[,-(1:4)],function(y){mean(abs((as.numeric(y)-u.OR)/s.OR),na.rm=T)})]
    )
    
    ### Calculate tMAD Score (slightly wrong)
    #No actual publication. Below is the only published material.
    #https://bgcs.org.uk/EMoore_BGCS_Oral.pdf
    
    x <-  cbind(detect$sWG$segments[,1:4,with=F],apply(detect$sWG$segments[,-(1:4)],2,function(y){y[y < -2] <- NA; y}))
    detect$t.MAD <- cbind(
      detect$sWG$meta,
      t.MAD=x[,sapply(.SD[,-(1:4)],function(y){mean(abs(y),na.rm=T)})]
    )
    
    saveRDS(detect,fname)
  }
  return(detect)
}

get.measures <- function(detect){
  return(cbind(
    detect$ichorCNA.optimal$summary[,.(ichorCNA.best=1-n_est)],
    detect$ichorCNA.low$summary[,.(ichorCNA.low=1-n_est)],
    detect$ichorCNA.default$summary[,.(ichorCNA.default=1-n_est)],
    detect$ichorCNA.suggested$summary[,.(ichorCNA.suggested=1-n_est)],
    detect$z.OR[,.(z.OR=z.OR)],
    detect$t.MAD[,.(t.Mad=t.MAD)],
    detect$absolute[,.(absolute=purity)]
  ))
}

get.measures.meta <- function(detect) return(cbind(detect$sWG$meta,get.measures(detect)))

detect.measures.output <- function(data,output.dir="~/Desktop/detect.data_dump/"){
  detect.output.all(data)
  fwrite(get.measures.meta(data),paste0(output.dir,"detect.measures.csv"))
}
