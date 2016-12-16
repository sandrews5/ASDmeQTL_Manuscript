#Here we make a volcano plot to make sure our meQTL discovery is working as intended. 
#We expect a degredation of significance with increasing distance. 

#Used to make Supplementary Figure 1. 

library(minfi)

resultsdirectory<-"" #Where the meQTL results are. 

chr<-22
myfiles<-list.files(path=resultsdirectory)
thischr<-grep(paste("restotal",chr,sep=""),myfiles,fixed=TRUE)
#Collect CpG site and SNP site position information and p-value. 
myres<-lapply(1:length(thischr),function(x){load(paste(resultsdirectory,myfiles[thischr[x]],sep=""))
			mypos<-pos[match(ret$gene,rownames(resids))]
			transP<-(-1*log(ret$pvalue,base=10))
			return(cbind(mypos,as.numeric(as.character(ret$snps)),transP))
})

#Format and plot. 
plotme<-data.frame(do.call("rbind",lapply(1:length(thischr),function(x){return(myres[[x]])})))
names(plotme)<-c("ProbePos","SNPPos","Pval")
plotme$Distance<-plotme$ProbePos-plotme$SNPPos
png("VolcanoPlot.png", width=1000, height=1000, bg="white", type="cairo") 
with(plotme, plot(Distance, Pval, pch=20, main=NULL, xlim=c(-1000000,1000000)))
dev.off()