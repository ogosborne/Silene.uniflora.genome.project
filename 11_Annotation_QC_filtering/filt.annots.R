library(wesanderson)
library(RColorBrewer)
#####################
# LOAD DATA
#####################

# names of different annotation stages: Evidence Modeler and three rounds of annotation refinement with PASA
DS.names <- c("1.all.evm","2.mydb_pasa.sqlite.gene_structures_post_PASA_updates.27063","3.mydb_pasa.sqlite.gene_structures_post_PASA_updates.34512","4.mydb_pasa.sqlite.gene_structures_post_PASA_updates.25334")

# load gffs and extract 5' + 3' info
gffs <- list()
gff.colnames <- c("seqid","source","type","start","end","score","strand","phase","attributes")
for(n in 1:4){
  # open gff
  gffs[[n]] <- read.delim(paste("data/gff3/",DS.names[[n]],".gff3",sep=""), header=F, comment.char="#",col.names = gff.colnames)
  # extract names of mRNAs with 3' and 5' UTRs
  utr5 <-  try(gsub("Parent=","",t(as.data.frame(strsplit(gffs[[n]][which(gffs[[n]]$type == "five_prime_UTR"),"attributes"],";")))[,2]),silent = T)
  utr3 <-  try(gsub("Parent=","",t(as.data.frame(strsplit(gffs[[n]][which(gffs[[n]]$type == "three_prime_UTR"),"attributes"],";")))[,2]),silent=T)
  # filter to keep only mRNAs
  gffs[[n]] <- gffs[[n]][which(gffs[[n]]$type == "mRNA"),c("seqid","source","start","end","strand","attributes")]
  gffs[[n]][which(gffs[[n]]$source == "."),"source"] <- "PASA"
  # get data from attribute line and split into gene, mRNA and name lines
  atts <- t(as.data.frame(strsplit(gffs[[n]]$attributes,";")))
  gffs[[n]]$mRNA <- gsub("ID=","",atts[,1])
  gffs[[n]]$gene <- gsub("Parent=","",atts[,2])
  gffs[[n]]$name <- gsub("%3B",";",
                         gsub("%3A",":",
                              gsub("%20"," ",
                                   gsub("Name=","",atts[,3]))))
  gffs[[n]]$attributes <- NULL
  rm(atts)
  # add utr presence
  if(class(utr5) != "try-error"){
    gffs[[n]]$five_prime_UTR <- gffs[[n]]$mRNA %in% utr5
  } else {
    gffs[[n]]$five_prime_UTR <- rep(FALSE,nrow(gffs[[n]]))
  }
  if(class(utr3) != "try-error"){
    gffs[[n]]$three_prime_UTR <- gffs[[n]]$mRNA %in% utr3
  } else {
    gffs[[n]]$three_prime_UTR <- rep(FALSE,nrow(gffs[[n]]))
  }
  rm(list=c("utr3","utr5"))
}
rm(list=c("n","gff.colnames"))

# load EVM evidence (made in script: evm.sh)
EVM.evidence <- read.csv("data/EVM_evidence/evidence.csv",header = T,stringsAsFactors = F)
EVM.evidence.RNA <- EVM.evidence[which(EVM.evidence$n.pasa > 0 | EVM.evidence$n.stringtie > 0),"gene.name"]
for(n in 1:4){
  gffs[[n]]$RNA_evidence <- FALSE
  gffs[[n]][which(gffs[[n]]$source == "PASA"),"RNA_evidence"] <- TRUE
  gffs[[n]][which(gffs[[n]]$gene %in% EVM.evidence.RNA),"RNA_evidence"] <- TRUE
}
rm(list=c("EVM.evidence.RNA","n"))

# load eggNOGmapper annotations (input files made in script: emapper.sh)
eggnog <- list()
eggnog.colnames <- c("mRNA","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
for(n in 1:4){
  eggnog[[n]] <- list()
  c = 70 # query coverage used
  # read emapper results
  eggnog[[n]][[paste("cov",c,sep="")]] <- read.table(paste("./data/eggNOGmapper/",DS.names[[n]],"_cov",c,"_e0.0001","/eggNOG_mapper.",DS.names[[n]],".emapper.annotations",sep=""),sep="\t",comment.char = "#",header = F,stringsAsFactors = F,na.strings = "-",quote = '"',col.names = eggnog.colnames)
  # just keep useful info
  eggnog.short <- eggnog[[n]][[paste("cov",c,sep="")]][,c("mRNA","evalue","Preferred_name","Description","COG_category","EC","KEGG_Pathway","PFAMs")]
  # add coverage to col names
  colnames(eggnog.short) <- c("mRNA",paste("eggnog_cov",c,"_",c("evalue","Preferred_name","Description","COG_category","EC","KEGG_Pathway","PFAMs"),sep=""))
  # add to output
  gffs[[n]] <- merge(gffs[[n]],eggnog.short,by="mRNA",all.x=T)
  rm(eggnog.short)
}
rm(list=c("n","c","eggnog.colnames"))

# load eggNOGmapper PFAM results (input files made in script: emapper.sh)
# transposon pfams
TE.pfams <- read.table("data/eggNOGmapper/transposon.pfams.tsv",sep="\t",header=F,stringsAsFactors = F,col.names = c("ID","ACC","Desc."))
pfam.colnames <- c("query_name","hit","evalue","sum_score","query_length","hmmfrom","hmmto","seqfrom","seqto","query_coverage")
# 
for(n in 1:4){
  c = 70 # query coverage used
  # read pfam results
  pfam <- read.table(paste("./data/eggNOGmapper/",DS.names[[n]],"_cov",c,"_e0.0001","/eggNOG_mapper.",DS.names[[n]],".emapper.pfam",sep=""),sep="\t",comment.char = "#",header = F,stringsAsFactors = F,na.strings = "",quote = '"',col.names = pfam.colnames)
    TE.genelist <- pfam[which(pfam$hit %in% TE.pfams$ID),"query_name"]
  # add to output
  gffs[[n]][,paste("eggnog_cov",c,"_TE_gene",sep="")] <- gffs[[n]]$mRNA %in% TE.genelist
    rm(list=c("TE.genelist","pfam"))
}
rm(list=c("n","c","pfam.colnames"))

# load lists of gene completeness for filtering (files created in script: GffReadFilt.sh)
# N: non-canonical splice site
# V: in-frame stop codons 
# U: single exon transcripts
# J: missing START or terminal STOP, or an in-frame stop codon
for(n in 1:4){
  for(f in c("J","N","U","V")){
    # load results
    gffread_filt <- readLines(paste("./data/gffread_filters/",DS.names[[n]],".gff3_IDs.pass.",f,".txt",sep=""))
    # add to output
    gffs[[n]][,paste("pass_filt_",f,sep="")] <- gffs[[n]]$mRNA %in% gffread_filt
    rm(gffread_filt)
  }
}
rm(list = c("n","f"))

# load busco results
busco.colnames <- c("Busco.id","Status","mRNA","Score","Length","OrthoDB url","Description")
busco <- list()
for(n in 1:4){
  busco[[n]] <- list()
  # for 3 different busco DBs
  for (ds in c("embryophyta","eudicots","viridiplantae")){
    busco[[n]][[ds]] <- read.table(paste("./data/BUSCO/",DS.names[[n]],"_",ds,"_odb10_full_table.tsv",sep=""),sep="\t",comment.char = "#",header = F,stringsAsFactors = F,na.strings = "",quote = '"',col.names = busco.colnames,fill = T)
    busco2add <- busco[[n]][[ds]][,c("mRNA","Busco.id","Status")]
    colnames(busco2add) <- c("mRNA",paste(ds,"Busco.id",sep="."),paste(ds,"Busco.status",sep="."))
    gffs[[n]] <-  merge(gffs[[n]],busco2add,by="mRNA",all.x=T)
    rm(busco2add)
  }
}
rm(list=c("n","ds","busco.colnames", "DS.names"))

#####################
# FUNCTONS FOR CHECKING AND FILTERING ANNOTATIONS
#####################

check_busco <- function(x,n.totals){
  # x: a data frame with at least the following columns: gene, <dataset>.Busco.status and <dataset>.Busco.id, where <dataset> corresponds to the name of each element in n.totals.
  # n.totals: a named vector with one element for each busco dataset in x, containing the total number of buscos in each dataset.
  #initialise output
  result <- list()
  # for each dataset
  for (n in 1:length(n.totals)){
    N.buscos <- rep(NA,4)
    names(N.buscos) <- c("Complete single-copy","Complete duplicated","Fragmented","Missing")
    # get all non-fragmented
    complete <- x[which(x[,paste(names(n.totals)[[n]],"Busco.status",sep=".")] == "Complete" | x[,paste(names(n.totals)[[n]],"Busco.status",sep=".")] == "Duplicated"),]
    # remove duplicates (different mRNAs from the same gene)
    complete <- complete[!duplicated(complete[,c("gene",paste(names(n.totals)[[n]],"Busco.id",sep="."))]),paste(names(n.totals)[[n]],"Busco.id",sep=".")]
    # frequency table of BUSCOs
    complete <- table(complete)
    # count complete single copy and duplicated
    N.buscos[["Complete single-copy"]] <- length(which(complete == 1))
    N.buscos[["Complete duplicated"]] <- length(which(complete > 1))
    rm(complete)
    # count fragmented
    N.buscos[["Fragmented"]] <- length(which(x[,paste(names(n.totals)[[n]],"Busco.status",sep=".")] == "Fragmented"))
    # count missing
    N.buscos[["Missing"]] <- n.totals[[n]]-sum(N.buscos[["Complete single-copy"]],N.buscos[["Complete duplicated"]],N.buscos[["Fragmented"]])
    # add to output
    result[[names(n.totals)[[n]]]] <- N.buscos
    rm(N.buscos)
  }
  result
}

# count complete genes
check_genes <- function(x){
  # x: a data frame with at least the following columns: gene, mRNA, pass_filt_J,  pass_filt_U, five_prime_UTR, three_prime_UTR. pass_filt_* columns are logicals indicating where the annotation passed the J, N, U and V filters from gffread, five_prime_UTR and three_prime_UTR are logicals indicating whether the annotation includes a 5' and 3' UTR.
  result <- rep(NA,8)
  names(result) <- c("N.genes","N.both.UTRs.mRNA.complete","N.3prime.UTR.mRNA.complete","N.5prime.UTR.mRNA.complete","N.no.UTR.mRNA.complete","N.mRNA.incomplete","N.multi.exon","N.single.exon")
  # gene lists
  genes <- list()
  # count genes
  result[["N.genes"]] <- length(unique(x$gene))
  # get list of genes
  genes[[1]] <- unique(x[which(x$pass_filt_J & x$three_prime_UTR & x$five_prime_UTR),"gene"])
  genes[[2]] <- unique(x[which(x$pass_filt_J & x$three_prime_UTR & x$five_prime_UTR == FALSE),"gene"])
  genes[[3]] <- unique(x[which(x$pass_filt_J & x$three_prime_UTR == FALSE & x$five_prime_UTR),"gene"])
  genes[[4]] <- unique(x[which(x$pass_filt_J & x$three_prime_UTR == FALSE & x$five_prime_UTR == FALSE),"gene"])
  genes[[5]] <- unique(x[which(x$pass_filt_J == FALSE),"gene"])
  genes[[6]] <- unique(x[which(x$pass_filt_U),"gene"])
  genes[[7]] <- unique(x[which(x$pass_filt_U == FALSE),"gene"])
  # count gene categories
  result[["N.both.UTRs.mRNA.complete"]] <- length(genes[[1]])
  result[["N.3prime.UTR.mRNA.complete"]] <- length(setdiff(genes[[2]],genes[[1]]))
  result[["N.5prime.UTR.mRNA.complete"]] <- length(setdiff(genes[[3]],c(genes[[1]],genes[[2]])))
  result[["N.no.UTR.mRNA.complete"]] <- length(setdiff(genes[[4]],c(genes[[1]],genes[[2]],genes[[3]])))
  result[["N.mRNA.incomplete"]] <- length(setdiff(genes[[5]],c(genes[[1]],genes[[2]],genes[[3]],genes[[4]])))
  result[["N.multi.exon"]] <- length(genes[[6]])
  result[["N.single.exon"]] <- length(setdiff(genes[[7]],genes[[6]]))
  rm(genes)
  result
}


#####################
# FILTER ANNOTATIONS
#####################
#
# FILTERS USED:
  #   - remove genes with non-canonical splice sites
  #   - remove genes with in-frame stop codons
  #   - removed genes missing START or terminal STOP
  #   - removed genes annotated as a TE gene by eggnogmapper
  #   - removed genes without EITHER RNA evidence OR a 70% coverage eggnogmapper hit

best <- gffs[[4]][which(gffs[[4]]$pass_filt_N & gffs[[4]]$pass_filt_V & gffs[[4]]$pass_filt_J & gffs[[4]]$eggnog_cov70_TE_gene == FALSE &  (gffs[[4]]$RNA_evidence | !is.na(gffs[[4]]$eggnog_cov70_evalue) )),]

# write mRNA IDs of best to file for and make gff with only those genes with gffread 
write(best$mRNA, file = "./results/best.IDs.txt")
system("gffread --ids best.IDs.txt --keep-genes 4.mydb_pasa.sqlite.gene_structures_post_PASA_updates.25334.gff3 > best.filt.set.gff3")
system("gffread -g Su_softmasked.fasta -x best.filt.annots.cds.fasta -y best.filt.annots.pep.fasta best.filt.annots.gff3")
# BUSCO plot
# get busco stats
busco_totals <- c(1614,2326,425)
names(busco_totals) <- c("embryophyta","eudicots","viridiplantae")
buscos_best <- do.call(cbind,check_busco(gffs[[4]],busco_totals ))
# convert to % and put in order
for(d in colnames(buscos_best)) buscos_best[,d] <- buscos_best[,d]/sum(buscos_best[,d])*100
buscos_best <- buscos_best[,c(2,1,3)]
colnames(buscos_best) <- c("Eudicots","Embryophyta","Viridiplantae")

# plot
pdf("./results/Busco_EMcov70_filtNV_noTE_RNA.pdf")
par(mar=c(5.1,6,4,0))
layout(matrix(c(1,1,2),ncol=3))
bus.pal <- wes_palette("Zissou1")[c(1,2,3,5)]
barplot(buscos_best,col=bus.pal,space = 0,ylab="% of BUSCOs",cex.names=1.5,cex.axis=1.5,cex.lab=1.5,yaxt="n")
axis(2, at=seq(0,100,10))
par(mar=c(5.1,0,4,0))
plot.new()
legend("left",fill=bus.pal,legend = rownames(buscos_best),cex=1.5,bty = "n")
dev.off()

## gene completeness plot
gene.counts <- do.call(rbind,
                       list(check_genes(gffs[[1]]), # all EVM genes
                            check_genes(gffs[[2]]), # PASA refinement 1
                            check_genes(gffs[[3]]), # PASA refinement 2
                            check_genes(gffs[[4]]), # PASA refinement 3
                            check_genes(gffs[[4]][which(gffs[[4]]$pass_filt_V & gffs[[4]]$pass_filt_N),]), # no in frame STOP or non-canonical splice-sites
                            check_genes(gffs[[4]][which(gffs[[4]]$pass_filt_V & gffs[[4]]$pass_filt_N & (gffs[[4]]$RNA_evidence | !is.na(gffs[[4]]$eggnog_cov70_evalue))),]), # RNA evidence or eggnog hit
                            check_genes(gffs[[4]][which(gffs[[4]]$pass_filt_V & gffs[[4]]$pass_filt_N & (gffs[[4]]$RNA_evidence | !is.na(gffs[[4]]$eggnog_cov70_evalue)) &  gffs[[4]]$eggnog_cov70_TE_gene == FALSE),]) # no transposon PFAMs
))

pdf("./results/Gene.completeness.pdf")
layout(matrix(c(1,1,1,2,2,
                3,3,3,4,4),byrow = T,ncol=5))
par(mar=c(10,15,4.1,0), mgp=c(4,1,0))
gene.pal <- c(rev(brewer.pal(4,name = "Greens")),"lightcoral")
barplot(t(gene.counts[,2:6]),space=0,names=c("EVM","PASA round 1","PASA round 2","PASA round 3","filter CDS","filter evidence","filter TE"),las=2,col=gene.pal,ylab="N genes")
plot.new()
par(mar=c(10,0,4.1,2.1), mgp=c(2,1,0))
legend("left",legend=c("Complete gene","5' UTR absent","3' UTR absent","Both UTR absent","Incomplete CDS"),fill = gene.pal,bty="n", cex = 1.5)
par(mar=c(10,15,4.1,0), mgp=c(4,1,0))
barplot(t(gene.counts[,7:8]),space=0,names=c("EVM","PASA round 1","PASA round 2","PASA round 3","filter CDS","filter evidence","filter TE"),las=2,col=wes_palette("Darjeeling2"),ylab="N genes")
plot.new()
par(mar=c(10,0,4.1,2.1), mgp=c(2,1,0))
legend("left",legend=c("Multi-exon genes","Single exon genes"),fill =wes_palette("Darjeeling2")[1:2],bty="n", cex = 1.5)

dev.off()

