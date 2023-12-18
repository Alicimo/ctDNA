library(data.table)
library(readxl)
source(paste0(project.dir,"src/load.sWG.R"))

detect.load <- function(project.dir){
  fname <- paste0(project.dir,"data/detect.raw.RDS")
  if(file.exists(fname)){
    data <- readRDS(fname)
  } else {
    data <- list()
    
    data$RECIST <- detect.load.RECIST.clinCalled(project.dir) # <---- called by Ramona
    #data$RECIST <- detect.load.tumour.sizes(project.dir) <---- calculated from volumes
  
    data$treatment <- detect.load.treatment(project.dir)
    
    data$sWG <- list()
    data$sWG$meta <- detect.load.sWG.meta(project.dir)
    data$sWG$cns <- detect.load.QDNA(project.dir)
    data$sWG$segments <- detect.load.DNAcopy(project.dir)
    if(nrow(data$sWG$meta[!(Sample.name %in% colnames(data$sWG$cns))])!=0) stop("Samples not found")
    
    i <- c(1:4,match(data$sWG$meta$Sample.name,colnames(data$sWG$cns)))
    data$sWG$cns <- data$sWG$cns[,i,with=F]
    data$sWG$segments <- data$sWG$segments[,i,with=F]
    
    data$genotype <- fread(paste0(project.dir,"data/meta/detect/detect.genotype.csv"))
    
    #data$NGTAS <- fread(paste0(project.dir,"data/R2b-Detect_all_mutations_master_file.txt"))
    #data$NGTAS.meta <- detect.generate.NGTAS.meta(data$NGTAS)
    #data <- detect.clean.NGTAS(data)
    
    data$clinical <- detect.load.clinical(project.dir)
    #data$clinical <- data$clinical[Patient.ID %in% data$sWG$meta$Patient.ID]
    
    if(0){
      x <- data.table(Patient.ID=1:323)
      x[,sWG:=Patient.ID %in% data$sWG$meta[Sample.type %in% c("ctDNA","BC")]$Patient.ID]
      x[,treatment:=Patient.ID %in% data$treatment$Patient.ID]
      x[,genotype:=Patient.ID %in% data$genotype$Patient.ID]
      x[,clinical:=Patient.ID %in% data$clinical$Patient.ID]
      x[,recist:=Patient.ID %in% data$RECIST$Patient.ID]
      x[x==FALSE] <- 0
      UpSetR::upset(x[,-1],order.by = "freq")
    }
    
    if(nrow(x <- data$genotype[data$clinical,on="Patient.ID"][(ER.status != pbc_er_status) | (Her2.status != pbc_her2_status)][,.(Patient.ID,ER.status,Her2.status,pbc_er_status,pbc_her2_status)])){
      warning(paste("Primary and metastatic genotype differs:",paste(x[,Patient.ID],collapse = ",")))
    }
    
    data <- detect.clean.sWG(data)
    
    pat.subset <- unique(data$sWG$meta[Sample.type %in% c("ctDNA","BC")]$Patient.ID)
    data$RECIST <- data$RECIST[Patient.ID %in% pat.subset]
    data$treatment <- data$treatment[Patient.ID %in% pat.subset]
    data$genotype <- data$genotype[Patient.ID %in% pat.subset]
    data$clinical <- data$clinical[Patient.ID %in% pat.subset]
    data <- subset.sWG(data,data$sWG$meta[,.I[Patient.ID %in% pat.subset]])
    
    saveRDS(data,fname)
  }
  return(data)
}  

detect.load.DNAcopy <- function(project.dir){
  dname <- paste0(project.dir,"data/sequencing/DNAcopy/")
  x <- lapply(list.files(dname,"*.tsv",full.names = T),fread)
  x <- do.call(cbind,c(x[[1]][,1:4],lapply(x,function(y) y[,5])))
  return(x)
}

detect.load.QDNA <- function(project.dir,bin.size=1000){
  dname <- paste0(project.dir,"data/sequencing/QDNA/bin.",bin.size,"/CNs/")
  x <- lapply(list.files(dname,"*.tsv",full.names = T),fread)
  x <- do.call(cbind,c(x[[1]][,1:4],lapply(x,function(y) y[,5])))
  return(x)
}

calculate.sWG.MAD <- function(x){
  y <- by(x,x$chromosome,function(y) apply(y[,-(1:4)], 2, diff))
  y <- do.call(rbind,y)
  y <- abs(y)
  y <- apply(y,2,median)
  return(y)
}

subset.sWG <- function(x,i){
  x$sWG$meta <- x$sWG$meta[i]
  x$sWG$cns <- x$sWG$cns[,c(1:4,i+4),with=F]
  x$sWG$segments <- x$sWG$segments[,c(1:4,i+4),with=F]
  print(paste("Sample count:",nrow(x$sWG$meta)))
  return(x)
}

detect.clean.sWG <- function(detect){

  print(paste("Sample count:",nrow(detect$sWG$meta)))
  
  #remove poor quality samples
  if(sum(i <- (calculate.sWG.MAD(detect$sWG$cns) > .2))){
    warning(paste("Poor quality samples:",paste(detect$sWG$meta[i,Sample.name],collapse = ",")))
    detect <- subset.sWG(detect,which(!i))
  }
  
  #remove samples with less than one million reads
  if(sum(i <- (detect$sWG$meta$aligned.reads < 1E6))){
    warning(paste("Samples lacking aligned reads:",paste(detect$sWG$meta[i,Sample.name],collapse = ",")))
    detect <- subset.sWG(detect,which(!i))
  }
  
  #remove unkwown samples
  if(sum(i <- detect$sWG$meta[,is.na(Sample.type)])){
    warning(paste("Samples with unknown types:",paste(detect$sWG$meta[i,Sample.name],collapse = ",")))
    detect <- subset.sWG(detect,which(!i))
  }
  
  #remove duplicate ctDNA with lower read counts
  i <- which( 
    (
      duplicated(detect$sWG$meta[,.(Sample.type,Patient.ID,Sample.num)]) |
      duplicated(detect$sWG$meta[,.(Sample.type,Patient.ID,Sample.num)],fromLast=T)
    ) &
      detect$sWG$meta[,Sample.type %in% c("ctDNA","BC")]
  )
  if(length(i)){
    warning(paste("Duplicated samples:",
                  paste(detect$sWG$meta[i][order(Patient.ID,Sample.type,Sample.num),Sample.name],
                        collapse = ",")))
    i <- i[ detect$sWG$meta[i][,.I[-(which.max(aligned.reads))],.(Patient.ID,Sample.type,Sample.num)][,V1] ]
    detect <- subset.sWG(detect,which(!(1:nrow(detect$sWG$meta) %in% i)))
  }
  
  return(detect)
}

detect.load.tumour.sizes <- function(project.dir){
  x <- fread(paste0(project.dir,"data/meta/detect/final.spreadsheet.csv"),check.names = T)
  y <- x[,1:7]
  z <- x[,-c(1:7,116:118)]
  
  q <- rbindlist(lapply(1:9,function(x){
    x <- z[,(1+(x-1)*12):((x*12)-1)]
    rbindlist(lapply(1:5,function(y){
      y <- melt(x,id.vars = c(1,2*y),measure.vars=(2*y)+1)
      names(y) <- c("date","location","scan.id","size")
      y
    }))
  }))
  q[q==""] <- NA
  q <- cbind(y[rep(1:dim(y)[1],dim(q)[1]/dim(y)[1]),],q)
  
  q[,scan.id:=tstrsplit(gsub("[i,v]","",scan.id),"\\.")[1]]
  q[,scan.id:=factor(scan.id,levels=c("B",paste0("X",1:8)),ordered = T)]
  q[,location:=toupper(location)]
  q[scan.id=="B",date:=closest.baseline]
  
  q <- q[!is.na(date)]
  q[location==0,location:=NA]
  q[is.na(location),size:=NA]

  q[,size:=as.numeric(size)]
  
  q <- q[order(patient.ID,scan.id)]
  q[!is.na(location),tumor.id:=1:.N,.(patient.ID,scan.id)]
  q[,lymph:=location %in% c("AX","HI","SC","EI")]
  
  q[lymph=="FALSE",target.les:=tumor.id %in% .SD[scan.id=="B" & size >= 10,tumor.id],patient.ID]
  q[lymph=="TRUE",target.les:=tumor.id %in% .SD[scan.id=="B" & size >= 15,tumor.id],patient.ID]

  
  #q[,scan.date.start:=as.Date(scan.date.start,"%d/%m/%Y")]
  #q[,scan.date.end:=as.Date(scan.date.end,"%d/%m/%Y") - scan.date.start]
  #q[,closest.baseline:=as.Date(closest.baseline,"%d/%m/%Y") - scan.date.start]
  
  q <- merge(q,detect.calc.response(q),sort = F)
  return(q)
}


detect.calc.response <- function(x){
  y <- x[,.(on.target.size=.SD[target.les=="TRUE",sum(size,na.rm = T)],
            on.target.size.nLN.CR=.SD[target.les=="TRUE" & lymph=="FALSE",sum(size,na.rm = T)==0],
            on.target.LN.CR=.SD[target.les=="TRUE" & lymph=="TRUE",all(size<10)],
            
            off.target.nLN.CR=.SD[target.les=="FALSE" & lymph=="FALSE",sum(size,na.rm = T)==0],
            off.target.LN.CR=.SD[target.les=="FALSE" & lymph=="TRUE",all(size<10)],
            
            off.target.tumours=list(.SD[target.les=="FALSE" & size!=0,tumor.id])),
         .(patient.ID,scan.id)]
  
  #print(y)
  
  y$response <- apply(y,1,function(z){ 
    p.ID <- z[1]
    s.ID <- z[2]
    onts <- as.numeric(z[3])
    onts.nln.CR <- z[4]
    onts.ln.CR <- z[5]
    offts.nln.CR <- z[6]
    offts.ln.CR <- z[7]
    offtt <- z[8]
    
    if (!("B" %in%  y[patient.ID==p.ID,scan.id])) return(NA)
    if (s.ID == "B") return(NA)
    
    if (y[patient.ID==p.ID & scan.id=="B",on.target.size==0]) onts.response <- NA
    else if (onts.nln.CR == "TRUE" & onts.ln.CR == "TRUE") onts.response <- "CR"
    else if ( (onts / y[patient.ID==p.ID & scan.id=="B",on.target.size]) <= .7 )  onts.response <- "PR"
    else if ( ((onts / y[patient.ID==p.ID & scan.id<s.ID,min(on.target.size)]) >= 1.2) 
              & (onts - y[patient.ID==p.ID & scan.id<s.ID,min(on.target.size)] >= 5 )) onts.response <- "PD"
    else onts.response <- "SD"
    
    if (offts.nln.CR == "TRUE" & offts.ln.CR == "TRUE") offts.response <- "CR"
    else if ( !(all(unlist(offtt) %in% unlist(y[patient.ID==p.ID & scan.id=="B",off.target.tumours]))))  offts.response <- "PD"
    else offts.response <- "Non-CR/Non-PD"
    
    #print(c(p.ID,s.ID,onts.response,offts.response))
    
    if(is.na(onts.response)) return(offts.response)
    if(onts.response == "CR" & offts.response == "CR") return("CR")
    if(onts.response == "PD" | offts.response == "PD") return("PD")
    if(onts.response %in% c("CR","PR") & offts.response != "PD") return("PR")
    if(onts.response == "SD" & offts.response != "PD") return("SD")
    else stop("Non-acceptable response categorisation")
  })
  
  y[,.(patient.ID,scan.id,response)]
}

detect.load.sWG.meta <- function(project.dir){
  x <- fread(paste0(project.dir,"data/meta/detect/sample.list.csv"))
  if(nrow(x[duplicated(Sample.name)])) stop("Duplicates with different DETECT status")
  
  x <- x[DETECT=="YES"]
  
  x[grepl("\\d+_{0,1}V\\d+",Sample.name),c("Patient.ID","Sample.num") := tstrsplit(regmatches(Sample.name,regexpr("\\d+_{0,1}V\\d+",Sample.name)),"V")]
  x[,Patient.ID:=gsub("_","",Patient.ID)]
  x[,Sample.num:=as.numeric(Sample.num)]
  x[!is.na(Patient.ID), Sample.type := "ctDNA"]
  x[is.na(Patient.ID), Patient.ID := regmatches(Sample.name,regexpr("\\d+",Sample.name))]
  x[,Patient.ID:=as.integer(Patient.ID)]
  x[grepl("BC",Sample.name), Sample.type := "BC"]
  x[grepl("PT",Sample.name), Sample.type := "PT"]
  x[grepl("MT",Sample.name), Sample.type := "MT"]
  
  y <- detect.load.ctDNA.dates(project.dir)
  x <- y[x,on=c("Patient.ID","Sample.num")]

  if(nrow(q <- x[Sample.type=="ctDNA" & is.na(Date.collected)])){
    stop(paste("Plasma samples with missing dates:",q[order(Patient.ID),paste(Patient.ID,Sample.num,sep=":",collapse = ",")]))
  }
  
  return(x)
}

detect.load.ctDNA.dates <- function(project.dir){
  fname <- paste0(project.dir,"data/meta/detect/detect.sample.log.xlsx")
  x <- data.table(read_xlsx(fname, sheet=1, skip=2,.name_repair = "minimal",trim_ws = T),check.names = T)
  x <- x[,.(Patient.ID=Study.number,
            Sample.num=as.integer(substr(Visit.number,2,3)),
            Date.collected=as.Date(Date.collected))]
  x <- unique(x)
  
  y <- fread(paste0(project.dir,"data/meta/detect/detect.extra_plasma_dates.csv"))
  y[,Date.collected:=as.Date(Date.collected,"%d/%m/%Y")]
  x <- rbind(x,y)
             
  if(nrow(q <- x[,.N,.(Patient.ID,Sample.num)][N>1])){
    stop(paste("Visits with multiple dates:",q[,paste(Patient.ID,Sample.num,sep=":",collapse = ",")]))
  }
  
  x <- x[order(Patient.ID,Sample.num)]
  
  if(nrow(i <- x[,.(diff(Date.collected),head(.I,-1)),Patient.ID][V1<=0])){
    stop(paste("Out of order dates:",x[i$V2,paste(Patient.ID,Sample.num,sep=":",collapse = ",")]))
    #warning(paste("Out of order dates:",x[i$V2,paste(Patient.ID,Sample.num,sep=":",collapse = ",")]))
    #x <- x[!(Patient.ID %in% i$Patient.ID)]
  }
  
  return(x)
}

detect.load.RECIST.clinCalled <- function(project.dir){
  x <- fread(paste0(project.dir,"data/meta/detect/detect.recist.csv"),check.names = T)
  if(nrow(w <- x[,.N,patient.ID][N>1])){
    #stop(paste("Multiple entries for patients:",w[,paste(patient.ID, collapse=",")]))
    warning(paste("Multiple entries for patients:",w[,paste(patient.ID, collapse=",")]))
    x <- x[,.SD[order(rowSums(.SD!=""))][1,],patient.ID]
  }
  
  y <- x[,1]
  z <- x[,-1]
  
  q <- rbindlist(lapply(1:11,function(i){
    x <- cbind(y,z[,(1+(i-1)*2):(i*2)])
    names(x) <- c("Patient.ID","date","RECIST")
    x$scan.id <- i
    return(x)
  }))
  q[q==""] <- NA
  q <- q[!is.na(date)]
  q[,date:=as.Date(date,"%d/%m/%Y")]
  q <- q[order(Patient.ID,date)]
  
  if(nrow(w <- q[,diff(scan.id),Patient.ID][V1!=1])){
    stop(paste("Patients with out of order dates:",w[,paste(unique(Patient.ID), collapse=",")]))
    warning(paste("Patients with out of order dates:",w[,paste(unique(Patient.ID), collapse=",")]))
    q <- q[!(Patient.ID %in% unique(w$Patient.ID))]
  }
  return(q)
}

detect.load.treatment <- function(project.dir){
  x <- fread(paste0(project.dir,"data/meta/detect/detect.treatment.csv"),check.names = T)
  y <- x[,1]
  z <- x[,-1]
  
  q <- rbindlist(lapply(1:(ncol(z)/4),function(i){
    x <- cbind(y,z[,(1+(i-1)*4):(i*4)])
    names(x) <- c("Patient.ID","treatment","start.date","end.date","notes")
    x$treatment.id <- paste0("X",i)
    return(x)
  }))
  q[q==""] <- NA
  q <- q[!(is.na(treatment) | is.na(start.date))]
  q <- q[order(Patient.ID,treatment.id)]
  
  #Two rows for one patient
  if(nrow(w <- q[,.N,.(Patient.ID,treatment.id)][N>1])){
    stop(paste("Duplicated patients:",w[,paste0(unique(Patient.ID),collapse = ",")]))
  }
  
  #Check start dates
  if(nrow(w <- q[!grepl("\\d{2}\\/\\d{2}\\/\\d{4}",start.date)])){
    stop(paste("Weird start dates found:",w[,paste0(Patient.ID,treatment.id,collapse = ",")]))
  }
  
  #Check end dates
  if(nrow(w <- q[!grepl("\\d{2}\\/\\d{2}\\/\\d{4}",end.date)][!is.na(end.date)])){
    warning(paste("Weird end dates found:",w[,paste0(Patient.ID,treatment.id,collapse = ",")]))
    q[!grepl("\\d{2}\\/\\d{2}\\/\\d{4}",end.date),end.date:=NA]
  }
  
  q[,start.date:=as.Date(start.date,"%d/%m/%Y")]
  q[,end.date:=as.Date(end.date,"%d/%m/%Y")]
  q[is.na(end.date),end.date:=Sys.Date()]
  q <- q[order(Patient.ID,start.date,end.date)]
  
  # End.date >= Start.date
  if(nrow(w <- q[end.date<start.date])){
    stop(paste("Weird end dates found:",w[,paste0(Patient.ID,treatment.id,collapse = ",")]))
  }
  
  # Weird treatment ordering
  if(nrow(w <- q[,.(treatment.id,c(diff(as.integer(substr(treatment.id,2,5))),NA)),Patient.ID][!is.na(V2)][V2!=1])){
    warning(paste("Weird treatment orders found:",w[,paste0(Patient.ID,treatment.id,collapse = ",")]))
  }
  
  # Single day treatments
  if(nrow(w <- q[start.date==end.date])){
    warning(paste("Single day treatments found:",w[,paste0(Patient.ID,treatment.id,collapse = ",")]))
    q <- q[!(start.date==end.date)]
  }
  
  # Merge overlapping treatments
  breaks <- q[, {
    tmp <- unique(sort(c(start.date, end.date)))
    .(start = head(tmp, -1L), end = tail(tmp, -1L))
  },Patient.ID]
  q <- q[breaks,
    on = .(Patient.ID, start.date <= start, end.date >= end),
   .(treatment=paste(sort(unique(treatment)),collapse=" + ")),
  by=.EACHI]
  q[treatment=="",treatment:="OFF"]

  # Remove "OFF periods if less than 7 days
  i <- q[,which(treatment=="OFF" & (end.date-start.date)<7)]
  q$end.date[i-1] <- q$start.date[i+1]
  q <- q[-i]
  
  return(q)
}

detect.load.clinical <- function(project.dir){
  x <- lapply(list.files(paste0(project.dir,"data/detect.clinical"),full.names = T),fread)
  x <- x[sapply(x,function(y) length(unique(y)[,.N,subject_id][,unique(N)]))==1]
  x <- Reduce(function(x, y) merge(x, unique(y),by="subject_id",all=T),x)
  x[,Patient.ID:=as.integer(substr(subject_id,3,10))]
  x[pbc_er_status=="Not Known",pbc_er_status:=NA]
  x[pbc_her2_status=="Not Known",pbc_her2_status:=NA]
  return(x)
  }

detect.clean.NGTAS <- function(detect){
  snames <- detect$NGTAS.meta[grepl("PD3",SampleID) & Sample.type=="ctDNA",SampleID]
  print("PD3 samples:")
  print(snames)
  detect$NGTAS.meta <- detect$NGTAS.meta[!(SampleID %in% snames)]
  detect$NGTAS <- detect$NGTAS[!(SampleID %in% snames)]
  
  dups <- detect$NGTAS.meta[Sample.type=="ctDNA",.(.N,list(SampleID)),.(Patient.ID,Sample.num)][N>=2]
  dups <- rbindlist(lapply(1:nrow(dups),function(i) detect$NGTAS[SampleID %in% dups[i,V2][[1]],mean(coverage),.(SampleID,Pool)][order(V1)][-.N] ))
  print("Lower quality dups:")
  print(dups$SampleID)
  detect$NGTAS.meta <- detect$NGTAS.meta[!(paste(SampleID,Pool) %in% dups[,paste(SampleID,Pool)])]
  detect$NGTAS <- detect$NGTAS[!(paste(SampleID,Pool) %in% dups[,paste(SampleID,Pool)])]
  
  detect$NGTAS <- rbindlist(
    lapply(unique(detect$NGTAS$patient),
           function(p.ID) detect$NGTAS[patient==p.ID & mutID %in% unique(detect$NGTAS[patient==p.ID & VAF >= .05, mutID])]
    )
  )
  
  return(detect)
}


detect.generate.NGTAS.meta <- function(dt.NGTAS){
  x <- unique(dt.NGTAS[,.(SampleID,Pool)])
  #if( nrow(x[duplicated(SampleID)]) != 0) stop("Duplicates") <----fix later
  
  x[grepl("\\d+_{0,1}V\\d+",SampleID),c("Patient.ID","Sample.num") :=
      tstrsplit(regmatches(SampleID,regexpr("\\d+_{0,1}V\\d+",SampleID)),"V")]
  x[grepl("\\d+_PLASMA_V\\d+",SampleID),c("Patient.ID","Sample.num") :=
      tstrsplit(regmatches(SampleID,regexpr("\\d+_PLASMA_V\\d+",SampleID)),"_PLASMA_V")]
  x[,Patient.ID:=gsub("_","",Patient.ID)]
  x[,Sample.num:=as.numeric(Sample.num)]
  x[!is.na(Patient.ID), Sample.type := "ctDNA"]
  x[is.na(Patient.ID), Patient.ID := regmatches(SampleID,regexpr("\\d+",SampleID))]
  x[,Patient.ID:=as.integer(Patient.ID)]
  x[grepl("BC",SampleID), Sample.type := "BC"]
  x[grepl("PT",SampleID), Sample.type := "PT"]
  x[grepl("MT",SampleID), Sample.type := "MT"]
  
  y <- detect.load.ctDNA.dates(project.dir)
  x <- merge(x,y,all.x=T)
  
  return(x)
}

detect.output.all <- function(data,output.dir="~/Desktop/detect.data_dump/"){
  dir.create(output.dir,showWarnings = F)
  fwrite(reshape(data$RECIST,idvar="Patient.ID",timevar="scan.id",direction="wide"),
         paste0(output.dir,"detect.recist.csv"))
  fwrite(reshape(data$treatment[,.(1:.N,start.date,end.date,treatment),Patient.ID]
                 ,idvar="Patient.ID",timevar="V1",direction="wide"),
         paste0(output.dir,"detect.treatment.csv"))
  fwrite(data$genotype,paste0(output.dir,"detect.genotype.csv"))
  fwrite(data$clinical,paste0(output.dir,"detect.NHS_open_clinica.csv"))
  fwrite(data$sWG$meta,paste0(output.dir,"detect.sWG_meta.csv"))
  fwrite(data$sWG$cns,paste0(output.dir,"detect.sWG_log2ratio.csv"))
  fwrite(data$sWG$segments,paste0(output.dir,"detect.sWG_segments.csv"))
} 
