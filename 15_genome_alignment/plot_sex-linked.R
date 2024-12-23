library(pafr)
library(stringr)
library(GenomicRanges)
# transcriptome mappings
mart.so <- read_paf("results/martin.So.paf")
mart.sp <- read_paf("results/martin.Sp.paf")
fila.sl <- read_paf("results/filatov.Sl.paf")
zemp.sl <- read_paf("results/zemp.Sl.paf")
# get sex-linked contig names
mart.info <- read.csv("./data/Martin_So_Sp_contig_info.csv",header=T)
mart.so.zw <- mart.info[which(mart.info$Species == "So" & mart.info$sexlinked == "ZW" ),"contig"]
mart.sp.xy <- mart.info[which(mart.info$Species == "Sp" & mart.info$sexlinked == "XY" ),"contig"]
fila.info <- read.csv("./data/Filatov_Sl_contig_info.csv",header=T)
fila.sl.xy <- fila.info[,"contig"]
zemp.info <- read.csv("./data/Zemp_Sl_contig_info.csv",header=T)
zemp.sl.xy <- zemp.info[,"contig"]
# add sex-linkage to mappings
mart.so$sex <- "autosomal"
mart.sp$sex <- "autosomal"
fila.sl$sex <- "autosomal"
zemp.sl$sex <- "autosomal"
mart.so[which(mart.so$qname %in% mart.so.zw),"sex"] <- "ZW"
mart.sp[which(mart.sp$qname %in% mart.sp.xy),"sex"] <- "XY"
fila.sl[which(fila.sl$qname %in% fila.sl.xy),"sex"] <- "XY"
zemp.sl[which(zemp.sl$qname %in% zemp.sl.xy),"sex"] <- "XY"
# chromosome sizes
chr.info <- data.frame(scaff.name=c("s_36","s_37","s_39","s_41","s_196","s_710","s_866","s_890","s_1327","s_12570","s_12571","s_12572"),length=c(73068074,40723422,76861226,54318462,47434262,32168091,72107204,61185580,51371477,53072158,51090366,68611950))
chr.info <- chr.info[order(chr.info$length,decreasing = T),]
row.names(chr.info) <- NULL
chr.info$chr.name <- paste("Chr",str_pad(1:12, 2, pad = "0"),sep="")
# change scaff names to chr names
for(n in 1:12){
  mart.so[which(mart.so$tname == chr.info[n,"scaff.name"]),"tname"] <- chr.info[n,"chr.name"]
  mart.sp[which(mart.sp$tname == chr.info[n,"scaff.name"]),"tname"] <- chr.info[n,"chr.name"]
  fila.sl[which(fila.sl$tname == chr.info[n,"scaff.name"]),"tname"] <- chr.info[n,"chr.name"]
  zemp.sl[which(zemp.sl$tname == chr.info[n,"scaff.name"]),"tname"] <- chr.info[n,"chr.name"]
}
# keep only mappings to chromosomes
mart.so.filt <- mart.so[which(mart.so$tname %in% paste("Chr",str_pad(1:12, 2, pad = "0"),sep="")),]
mart.sp.filt <- mart.sp[which(mart.sp$tname %in% paste("Chr",str_pad(1:12, 2, pad = "0"),sep="")),]
fila.sl.filt <- fila.sl[which(fila.sl$tname %in% paste("Chr",str_pad(1:12, 2, pad = "0"),sep="")),]
zemp.sl.filt <- zemp.sl[which(zemp.sl$tname %in% paste("Chr",str_pad(1:12, 2, pad = "0"),sep="")),]
# add midpoint of mapping
mart.so.filt$midpoint <- (mart.so.filt$tstart+mart.so.filt$tend)/2
mart.sp.filt$midpoint <- (mart.sp.filt$tstart+mart.sp.filt$tend)/2
fila.sl.filt$midpoint <- (fila.sl.filt$tstart+fila.sl.filt$tend)/2
zemp.sl.filt$midpoint <- (zemp.sl.filt$tstart+zemp.sl.filt$tend)/2
# add numeric chr names
mart.so.filt$chr.num <- as.numeric(substr(mart.so.filt$tname,4,6))
mart.sp.filt$chr.num <- as.numeric(substr(mart.sp.filt$tname,4,6))
fila.sl.filt$chr.num <- as.numeric(substr(fila.sl.filt$tname,4,6))
zemp.sl.filt$chr.num <- as.numeric(substr(zemp.sl.filt$tname,4,6))
# add query coverage
mart.so.filt$qcov <- mart.so.filt$alen/mart.so.filt$qlen 
mart.sp.filt$qcov <- mart.sp.filt$alen/mart.sp.filt$qlen 
fila.sl.filt$qcov <- fila.sl.filt$alen/fila.sl.filt$qlen 
zemp.sl.filt$qcov <- zemp.sl.filt$alen/zemp.sl.filt$qlen 
# keep only mappings with query coverage >= 0.8 and mapping quality > 25 
mart.so.filt <- mart.so.filt[which(mart.so.filt$qcov >= 0.8 & mart.so.filt$mapq > 25),]
mart.sp.filt <- mart.sp.filt[which(mart.sp.filt$qcov >= 0.8 & mart.sp.filt$mapq > 25),]
fila.sl.filt <- fila.sl.filt[which(fila.sl.filt$qcov >= 0.8 & fila.sl.filt$mapq > 25),]
zemp.sl.filt <- zemp.sl.filt[which(zemp.sl.filt$qcov >= 0.8 & zemp.sl.filt$mapq > 25),]
zemp.sl.filt$zd <- NULL
# add to list
my.filt <- list(sl = rbind(fila.sl.filt,zemp.sl.filt),
                so =  mart.so.filt,
                sp = mart.sp.filt)
# function to plot sex-linked contigs position on S. uniflora chromosomes
plotSexLinked <- function(x,sex="XY",pos.col="midpoint",sex.col="sex",n.chrom=12,chrom.num.col="chr.num",aut.rgb=c(0.8,0.8,0.8),sex.colour="red",main="",xlab="Position (Mb)"){
  plot(c(0,max(x$tlen)/1000000),c(1,n.chrom),col="white",xlab=xlab,ylab="Chromosome",yaxt="n",main=main,bty='n')
  for (n in 1:nrow(x)){
    lines(rep(x[n,pos.col]/1000000,2),c(x[n,chrom.num.col]-0.4,x[n,chrom.num.col]+0.4),col=rgb(aut.rgb[[1]],aut.rgb[[2]],aut.rgb[[3]],0.1))
  }
  for (n in 1:nrow(x)){
    if(x[n,sex.col] == sex){
      lines(rep(x[n,pos.col]/1000000,2),c(x[n,chrom.num.col]-0.4,x[n,chrom.num.col]+0.4),col=sex.colour)
    }
  }
  mtext(1:n.chrom,2,at=1:n.chrom,las=2,line=1,cex=0.8)
  legend("topright",legend=c(sex,"autosomal","none"),fill=c(sex.colour,rgb(aut.rgb[[1]],aut.rgb[[2]],aut.rgb[[3]],1),"white"),bty="n")
}
# plot
pdf("./results/all.spp.trans.pdf")
par(mfrow=c(3,1),mar=c(3,4,3,2))
plotSexLinked(mart.so.filt,sex="ZW",main="S. otites",xlab="")
plotSexLinked(mart.sp.filt,sex="XY",main="S. pseudotites",xlab="")
plotSexLinked(both.sl.filt,sex="XY",main="S. latifolia",xlab="Position (Mb)")
dev.off()

