x <- detect$treatment[start.date!=end.date]

y <- rbindlist(lapply(1:nrow(x),function(i){
  tb <- x[i]
  
  subset.ctDNA <- merge(
    detect$sWG.meta[Sample.type=="ctDNA"][Patient.ID==tb$patient.ID & Date.collected>tb$start.date & as.Date(Date.collected)<tb$end.date][order(as.Date(Date.collected))],
    detect$ichorCNA.summary.default
  )
  
  subset.CT <- detect$RECIST$EMMA[patient.ID==tb$patient.ID & as.Date(date,"%d/%m/%Y") > tb$start.date & as.Date(date,"%d/%m/%Y") < tb$end.date][order(as.Date(date,"%d/%m/%Y"))]
  subset.CT[,date:=as.Date(date,"%d/%m/%Y")]
  
  subset.ctDNA[,NN:=which.min(abs(Date.collected-subset.CT$date)),Sample.name]
  subset.CT[,NN:=which.min(abs(subset.ctDNA$Date.collected-date)),scan.id]
  
  subset.ctDNA <- subset.ctDNA[subset.ctDNA[,subset.CT[NN]$NN==.I]]
  subset.ctDNA <- cbind(subset.ctDNA,subset.CT[subset.ctDNA$NN])
  
  subset.ctDNA[,.(Patient.ID,ctDNA.date=Date.collected,CT.date=date)]
}))
