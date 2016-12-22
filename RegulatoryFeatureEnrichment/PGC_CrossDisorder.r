#In this script, we will download and format the PGC cross-disorder results for 
#further use. 

#These results can be downloaded from here:

#Set directory and read in results.
setwd("/CrossDisorder")
f<-read.table("pgc.cross.full.2013-03.txt",header=TRUE,stringsAsFactors=FALSE,na.strings=".")

#Create a file to put into liftOver to get hg19 positions
hg18pos<-matrix(0,nrow(f),1)
hg18pos[,1]<-paste0("chr",f$hg18chr,":",f$bp,"-",f$bp)
write.table(hg18pos,file="/dcs01/feinberglab/roadmap/meQTL_Paper/CrossTissue/CrossDisorder/hg18pos.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Read in results output from liftOver software and format. 
library(rtracklayer)
hg19pos<-read.table("hglft_genome_1ee8_47b450.bed",stringsAsFactors=FALSE)
hg19pos<-unlist(lapply(hg19pos[,1],function(x){strsplit(x,":")[[1]][2]}))
hg19pos<-unlist(lapply(hg19pos,function(x){strsplit(x,"-")[[1]][1]}))

#Remove those positions that failed to map over to new genome coordinates.
#Label appropriate columns and name them in similar way to PGC ASD result object
#used earlier so a lot of the code can be adopted more readily. 
crosspgc<-f[-which(is.na(f$pos)),]
crosspgc<-crosspgc[,c("snpid","hg18chr","pos","or","pval","CEUaf")]
colnames(crosspgc)<-c("SNP","CHR","BP","OR","P","MAF")

#Save these results. 
save(crosspgc,file="/CrossDisorder/crosspgc.rda")
