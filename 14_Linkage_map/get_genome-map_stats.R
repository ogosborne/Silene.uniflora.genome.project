library(tidyverse)
##### Read data
# read map
mapdir <- "Linkage_map/"
map <- read.table(paste0(mapdir, "map.tsv"), header = T)
# read fasta index for sequence lengths
chr.info <- read.table("data/Su_softmasked.fasta.fai")[,1:2]
colnames(chr.info) <- c("scaff.name", "length")
# order by size
chr.info <- chr.info[order(chr.info$length, decreasing = T),]
rownames(chr.info) <- NULL
# rename top 12 scaffolds as Chr1 - 12
chr.info$chr.name <- chr.info$scaff.name
chr.info[1:12, "chr.name"] <-  paste("Chr",str_pad(1:12, 2, pad = "0"),sep="")

##### Assembly stats

# total assembly size
total_len <- sum(chr.info$length)
total_len
## 1269363085
cat("Total assembly length =", round(total_len / 1E6, 2), "Mb\n")
## Total assembly length = 1269.36 Mb
cat("Total assembly length =",round(total_len / 1E9, 2), "Gb\n")
## Total assembly length = 1.27 Gb

# N scaffolds
Nscaff <- nrow(chr.info)
Nscaff

# L50
L50 <- min(which(cumsum(chr.info$length) > total_len*0.5))
cat("L50 =", L50, "\n")
# L50 = 11

# N50
N50 <- chr.info[L50, "length"]
N50
## 40723422
cat("N50 =", round(N50 / 1E6, 2), "Mb\n")
## N50 = 40.72 Mb

# size of first 12 scaffs
top12_len <- sum(chr.info[1:12, "length"])
top12_len/1E6
# 682.0123
cat(round(top12_len/total_len*100, 2), "% of the assembly is contained in the largest 12 scaffolds")
## 53.73 % of the assembly is contained in the largest 12 scaffolds

# size range of top 12 scaffolds
round(range(chr.info[1:12, "length"])/1E6, 2)

# N scaffolds over 1Mb
Nover1Mb <- length(which(chr.info$length > 1E6))
Nover1Mb - 12

# N scaffolds in map
NscaffMap <- length(unique(map$scaffold))
NscaffMap

# scaffold-LG combination counts
LGscaff <- paste0(map$LG, map$scaffold) %>%
  table() %>%
  sort(decreasing=T)
LGscaff
# N markers
Nmar <- sum(LGscaff)
Nmar

# proportion of markers which match top LG
chr.info$LG <- NA
chr.info$N.mar <- 0
chr.info$N.mar.correct.LG <- 0
chr.info$Prop.mar.correct.LG <- NA
for(i in 1:nrow(chr.info)){
  tmp <- list()
  tmp$sc <- chr.info[i,"scaff.name"]
  tmp$ch <- chr.info[i,"chr.name"]
  tmp$df <- map[which(map$scaffold == tmp$sc), ]
  tmp$nm <- nrow(tmp$df)
  tmp$nt <- sort(table(tmp$df$LG), decreasing = T)[1]
  tmp$lg <- names(tmp$nt)
  tmp$pc <- ifelse(tmp$nm > 0, tmp$nt/tmp$nm, NA)
  if(tmp$nm > 0){
    chr.info[i,"LG"] <- tmp$lg
    chr.info[i,"N.mar"] <- tmp$nm
    chr.info[i,"N.mar.correct.LG"] <- tmp$nt
    chr.info[i,"Prop.mar.correct.LG"] <- tmp$pc
  }
  rm(tmp)
}

# mean percent of chrom markers which map to correct LG
round(mean( head(chr.info, 12)$Prop.mar.correct.LG)*100,2)
# 95.6
# range of percent of chrom markers which map to correct LG
round(range( head(chr.info, 12)$Prop.mar.correct.LG)*100,2)
# 87.20 97.45
