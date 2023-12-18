
### Find genes within amplifications/deletions

options(stringsAsFactors=F)

# working directory
setwd("/lustre/cclab/cliffo01/4Meiling/genelist")

# qdnaseq segments file
qdnaseq <- read.delim("24-F.duprem.bam_QDNAseq.txt")

# biomart hg19/GRCh37 with chrom, gene start/end, common name, ensembl id
mart <- read.csv("GRCh37-hg19_mart_export.txt")

# threshold for calling a amp/del
threshold <- 0.5


##########


# subset qdnaseq by threshold
segments <- qdnaseq[ , grep( "segments" , colnames(qdnaseq) ) ]
sub_qdnaseq <- qdnaseq[ which( abs(segments) > threshold ) , ]

# loop through subset qdnaseq identifying genes
genelist <- c()
for(i in 1:nrow(sub_qdnaseq)){

      # extract locus information
      tmprow <- sub_qdnaseq[i,]
      tmpchr <- as.character(tmprow$chromosome)
      tmpstart <- as.numeric(tmprow$start)
      tmpend <- as.numeric(tmprow$end)

      # which rows of biomart overlap with this
      mart_idx <- which(
      	       	  # same chromosome
      		  mart$Chromosome.Name == tmpchr &
		  ((
		    # the end or start of the segment lies within gene
		    (mart$Gene.Start..bp. < tmpend &
		  	mart$Gene.End..bp. > tmpend) |
		    (mart$Gene.End..bp. > tmpstart &
		    	mart$Gene.Start..bp. < tmpstart)
		  )|(
		    # or the segment encompasses the entire gene
		    mart$Gene.Start..bp. > tmpstart &
		    	mart$Gene.End..bp. < tmpend
		  ))
		  )

      # add rows to final genelist
      genelist <- rbind(genelist,mart[mart_idx,])

}

# write results to file
write.table(genelist,"genelist_output.txt",
		quote=F,sep="\t",row.names=F,col.names=T)


