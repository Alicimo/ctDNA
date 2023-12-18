root.dir <- "~/OneDrive/projects/"
project.dir <- paste0(root.dir,"ctDNA/")
source(paste0(project.dir,"src/load.detect.R"))

detect.flist <- detect.load.sWG.meta(project.dir)
sWG.flist <- colnames(sWG.load.data(paste0(project.dir,"data/sWG/1Mb/"))$cns)[-(1:4)]
i <- which(sWG.flist %in% detect.flist$Sample.name)

d <- paste0(project.dir,"data/sWG/1Mb/")
for(j in i){
  CMD <- paste("pdfseparate",
               "-f",j,
               "-l",j,
               paste0(d,"QDNA.merged.plot.pdf"),
               paste0(d,"indv.detect.plots/",sWG.flist[j],".pdf")
               )
  system(CMD)
}
