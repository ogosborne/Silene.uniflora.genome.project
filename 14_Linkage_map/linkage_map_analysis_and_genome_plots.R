library(stringr)
library(tidyverse)
library(dplyr)
library(GGally)
library(ggbio)
library(GenomicRanges)

##### 1. LOAD DATA #####
# load map
map <- read.table("Linkage_map/map.tsv", header = T)
# read fasta index for sequence lengths
chr.info <- read.table("data/Su_softmasked.fasta.fai")[,1:2]
colnames(chr.info) <- c("scaff.name", "length")
# order by size
scaff.info <- scaff.info[order(scaff.info$length, decreasing = T),]
rownames(scaff.info) <- NULL
# rename top 12 scaffolds as Chr1 - 12
scaff.info$chr.name <- scaff.info$scaff.name
scaff.info[1:12, "chr.name"] <-  paste("Chr",str_pad(1:12, 2, pad = "0"),sep="")
# extract chromosomes
chr.info <- scaff.info[1:12,]

##### 2. FILTER MARKERS WITH ASSEMBLY - MAP INCONSISTENCIES #####
# count scaffold-LG associations, to identify correct LG for each scaffold
scaffLGcounts <- map %>%
  count(scaffold, LG) %>%
  complete(scaffold, nesting(LG), fill = list(n = 0)) %>% 
  as.data.frame()
# reshape to counts table
scaffLGcounts <-  reshape(scaffLGcounts, idvar = "scaffold", timevar = "LG", direction = "w")  
colnames(scaffLGcounts) <- gsub("n.", "", colnames(scaffLGcounts), fixed = T)
# get main LG per scaffold
scaffLGcounts$mainLG <- colnames(scaffLGcounts)[2:ncol(scaffLGcounts)][max.col(scaffLGcounts[,2:ncol(scaffLGcounts)], ties.method = "f")]
# get N markers per scaffold
scaffLGcounts$nMarkers <- rowSums(scaffLGcounts[,2:(ncol(scaffLGcounts)-1)])
# add chr name and scaffold length
scaffLGcounts$Chr <- scaff.info[match(scaffLGcounts$scaffold, scaff.info$scaff.name), "chr.name"]
scaffLGcounts$length <- scaff.info[match(scaffLGcounts$scaffold, scaff.info$scaff.name), "length"]
# sort by chr length
scaffLGcounts <- scaffLGcounts[order(scaffLGcounts$length, decreasing = T),]
rownames(scaffLGcounts) <- scaffLGcounts$scaffold
# add Chr names
map$Chr <- scaff.info[match(map$scaffold, scaff.info$scaff.name), "chr.name"]
# mark markers on the wrong LG for the scaffold
map$LGokay <- map$LG == scaffLGcounts[map$scaffold,"mainLG"]
# get proportions of markers with correct LG per scaff
LGOkay_cnt <- table(map[,c("Chr", "LGokay")]) %>% 
  as.data.frame.matrix() 
LGOkay_prp <- LGOkay_cnt %>% 
  as.matrix() %>%
  proportions(margin = 1)
scaffLGcounts$prop_main_LG <- LGOkay_prp[scaffLGcounts$Chr,"TRUE"]
write.csv(scaffLGcounts, file = "results/scaffLGcounts.csv", row.names = F)
# identify blocks of consecutive markers from the same scaffold for filtering windows in recombination rate estimation.
blocks <- map %>% 
  group_by(id = consecutive_id(scaffold, LG), LG, scaffold, Chr, LGokay) 
blocks <- blocks[which(blocks$LGokay == TRUE & grepl(pattern = "Chr", blocks$Chr)), c("LG","scaffold", "Chr", "pos", "map.pos", "id")] %>%
  as.data.frame()
# remove outliers from map using loess regression
remove_map_outliers <- function(x, chr.info, loess.span = 0.1, ol.lim = 1.5){ 
  out <- list()
  for(i in 1:nrow(chr.info)){
    # get chr data
    scaff_name <- chr.info[i,"scaff.name"]
    chr_name <- chr.info[i,"chr.name"]
    chr_len_mb <- chr.info[i,"length"]/1E6
    my.df <- x[which(x$scaffold == scaff_name),]
    my.df$pos.mb <- my.df$pos/1E6
    # change orientation so cor is +ve
    my.cor <- cor(my.df$pos.mb, my.df$map.pos, method = "s")
    if(my.cor < 0) my.df$map.pos <- max(my.df$map.pos) - my.df$map.pos
    # filter marker for outliers using loess
    lo <- loess(my.df$map.pos ~ my.df$pos.mb, span = loess.span, control = loess.control(surface = "direct"))  
    my.df$loess.diff <- abs(my.df$map.pos - predict(lo, my.df$pos.mb))
    Q <- quantile(my.df$loess.diff, probs=0.75, na.rm = F) + IQR(my.df$loess.diff) * ol.lim
    my.df$outlier <- TRUE
    my.df[which(my.df$loess.diff < Q), "outlier"] <- FALSE
    my.df <- my.df[my.df$outlier == FALSE, c("LG", "scaffold", "Chr", "pos", "pos.mb", "map.pos", "id")]
    out[[chr_name]] <- my.df
  }
  out
}
blocks_noOL <- remove_map_outliers(blocks, chr.info = chr.info)

##### 3. CALCULATE RECOMBINATION RATE
# function to get recombination rate estimates
get.rr <- function(x, chr.info, wind.size = 1, min.nmarkers = 5, loess.span = 0.1,  max.rr = Inf){ 
  # set up output
  out <- list()
  for(i in names(x)){
    # get chr data
    chr_len_mb <- chr.info[which(chr.info$chr.name == i),"length"]/1E6
    my.df <- x[[i]]
    nmarkers <- nrow(my.df)
    # loess regression
    lo <- loess(my.df$map.pos ~ my.df$pos.mb, span = loess.span, control = loess.control(surface = "direct"))
    # get windows
    my.out <- data.frame(Chr = i, start = seq(from = 0, to = chr_len_mb, by = wind.size),
                         end = c(seq(from = wind.size, to = chr_len_mb, by = wind.size), chr_len_mb))
    my.out$mid <- (my.out$start + my.out$end) / 2
    my.out$rr.lo <- NA
    # get slope for windows
    for(w in 1:nrow(my.out)){
      my.start <- my.out[w,"start"]
      my.end <- my.out[w,"end"]
      my.wind.df <- my.df[which(my.df$pos.mb > my.start & my.df$pos.mb <= my.end),]
      # set window to NA if it contains markers from multiple scaffolds(i.e. contains more than one block) or has too few markers
      if(nrow(my.wind.df) < min.nmarkers | length(unique(my.wind.df$id)) > 1){
          my.rr.lo <- NA
      # otherwise get slope
      } else {
          my.rr.lo <- abs(diff(predict(lo, c(my.start, my.end)))) / (my.end-my.start)
      }
      # set window to NA if it has recomb rate > max.rr
      if(!is.na(my.rr.lo)){
        if(my.rr.lo > max.rr) my.rr.lo <- NA
      } 
      my.out[w,"rr.lo"] <- my.rr.lo
    }
    out[[i]] <- my.out
  }
  out
}
# run
recomb <- get.rr(x = blocks_noOL, chr.info = chr.info, wind.size = 1, min.nmarkers = 1, loess.span = 0.1, max.rr = Inf)
# rr boxplot
recombs <- list()
for(i in names(recomb)){
  recombs[[i]] <- recomb[[i]]$rr.lo
}
chrom.cols = c(Chr01="firebrick",
               Chr02="red",
               Chr03="salmon",
               Chr04="darkorange",
               Chr05="gold",
               Chr06="yellowgreen",
               Chr07="forestgreen",
               Chr08="#449293",
               Chr09="steelblue",
               Chr10="cornflowerblue",
               Chr11="slateblue",
               Chr12="#571B90")
boxplot(recombs, col = chrom.cols, ylab = "Recombination rate (cM/Mb)")
# mean recomb. rate
chr.info$mean.rr <- lapply(recombs, mean, na.rm = T) %>% unlist()
chr.info$length.mb <- chr.info$length/1E6
cor.test(chr.info$length.mb, chr.info$mean.rr, method = "s")
plot(chr.info$length.mb, chr.info$mean.rr, col = chrom.cols, pch = 19, cex = 1.5, xlab = "Chromosome length (Mb)", ylab = "Recombination rate (cM/Mb)")
abline(lm(mean.rr ~ length.mb, data = chr.info))
# plot marey map with recombination rate on alt axis
pdf("results/marey_maps.pdf")
par(mfrow=c(4,3), mar = c(4,4,1,2), mgp = c(2,1,0))
for(i in names(recomb)) {
  my.df <- blocks_noOL[[i]]
  plot(my.df$pos.mb, my.df$map.pos, xlab = "Physical position (Mb)", ylab = "Map position (cM)", main = i, ylim = c(0,max(my.df$map.pos, na.rm = T)))
  #points(my.df.ol$pos.mb, my.df.ol$map.pos, col = "white")
  my.df <- recomb[[i]]
  my.df <- my.df[!is.na(my.df$rr.lo),]
  points(my.df$mid, my.df$rr.lo*5, col = "cornflowerblue", type = "p", pch = 19)
  axis(labels = seq(from = 0, to = 16, by = 2), at = c(0,seq(from=10, to = 80, by = 10)), side = 4,  las = 1, cex = 0.5, col = "cornflowerblue", col.axis = "cornflowerblue")
  mtext(text = "r (cM/Mb)", side = 4, line =-1, col = "cornflowerblue", cex = 0.7, las = 3)
}
dev.off()
# save
save(recomb,file = "results/recombination.rate.RData")

##### 4. PLOT LARGE SCAFFOLDS ON LGS
# scaffs over 1Mb which aren't chromosomes
# s_38;    13,055,758 bp; LG.1; Chr11
# s_12569;  3,394,578 bp; LG.4; Chr12
# s_7457;   2,760,734 bp; LG.3; Chr02
# s_1303;   1,724,594 bp; LG.4; Chr12
# s_10375;  1,072,322 bp; LG.4; Chr12
# plot LGs
# function to plot LGs
plot_scaff_pos <- function(map, scaffs, LG, scaff.info, cols = c("black", "royalblue", "lightcoral", "darkgoldenrod")){
  # get data for each scaffold
  dat <- list()
  for(s in scaffs){
    my.df <- map[which(map$LG == LG & map$Chr == s),]
    # change orientation so cor is +ve
    my.cor <- cor(my.df$pos, my.df$map.pos, method = "s")
    if(my.cor < 0){
      print(paste0("changing orientation: ", s))
      my.df$pos <- max(my.df$pos) - my.df$pos
    } 
    # order
    my.df <- my.df[order(my.df$pos),]
    dat[[s]] <- my.df
  }
  # get LG-wide rr
  my.LG <- map[which(map$LG == LG & map$LGokay == TRUE),]
  LGscaffs <- unique(my.LG$scaffold)
  LG_physlen <- sum(scaff.info[scaff.info$scaff.name %in% LGscaffs, "length"])
  LG_maplen <- max(my.LG$map.pos)
  rr <- LG_maplen/LG_physlen
  # add buffer to scaffold positions
  sc.pos <- c()
  for(s in c(scaffs)){
    sc.pos[s] <- dat[[s]][1, "map.pos"]
  }
  sc.pos <- sc.pos / rr
  for(s in c(scaffs)){
    dat[[s]]$pos <- dat[[s]]$pos + sc.pos[s]
  }
  # name colours
  cols <- cols[1:length(scaffs)]
  names(cols) <- scaffs
  # combine data from each scaffold
  dat <- Reduce(rbind, dat)
  # plot
  plot(dat$pos/1E6, dat$map.pos, 
       col = cols[dat$Chr], 
       xlab = "Physical position (Mb)", 
       ylab = "Map position (cM)", 
       main = LG,
       pch = 19)
  legend("bottomright", legend = names(cols), col = cols, pch = 19, bty = "n")
}
# plot
plot_scaff_pos(map = map, scaffs = c("Chr11", "s_38"), LG = "LG.1", scaff.info = scaff.info)
plot_scaff_pos(map = map, scaffs = c("Chr02", "s_7457"), LG = "LG.3", scaff.info = scaff.info)
plot_scaff_pos(map = map, scaffs = c("Chr12", "s_12569", "s_1303", "s_10375"), LG = "LG.4", scaff.info = scaff.info)

#### 5. GENOMIC FEATURE RELATIONSHIPS
# load data for plot
scaff2chrom <- chr.info$chr.name
names(scaff2chrom) <- chr.info$scaff.name
chrom2scaff <- chr.info$scaff.name
names(chrom2scaff) <- chr.info$chr.name
# Gene density
NGenes <- read.table("../../../results/14_genome_structure_stats/gene_density/ngenes_1Mb.tsv")
colnames(NGenes) <- c("scaffold", "start", "end", "N.genes")
# GC content
GC <- read.table("../../../results/14_genome_structure_stats/gc/GC_1Mb.tsv")
colnames(GC) <- c("scaffold", "start", "end", "pct_at", "pct_gc",	"num_A", "num_C", "num_G", "num_T", "num_N", "num_oth", "seq_len")
GC <- GC[,c("scaffold", "start", "end", "pct_gc")]
# Repeat content
reps <- read.table("data/RepeatMasker/repeats.gff", sep = "\t")[,c(1,4,5,9)] # repeats.gff created in RepeatMasker.sh
colnames(reps) <- c("scaff", "start", "end", "info")
reps$scaff <- paste0("s_", reps$scaff)
reps <- reps[which(reps$scaff %in% chr.info$scaff.name),]
# split seq names to get 
inf <- str_split_fixed(reps$info, pattern = "\\|", n = 7)
# get repeat family
repfam <- ifelse(inf[,3] == "", "simple", inf[,3])
# simplify repeat family classifications
rf_c <- c(rep("DNA transposon", 10),
          "non-LTR retrotransposon",
          "Other LTR retrotransposon",
          "Copia",
          "Gypsy",
          "Other/Unknown",
          "non-LTR retrotransposon",
          "Other/Unknown",
          "Simple repeat",
          "DNA transposon",
          "Other/Unknown",
          "rRNA",
          "Satellite",
          "Simple repeat",
          "DNA transposon")
names(rf_c) <- c("DNA","DNA/En-Spm","DNA/Harbinger","DNA/hAT","DNA/hAT-Ac","DNA/Mite","DNA/MuDR","DNA/Stowaway","DNA/TcMar","DNA/Tourist","LINE","LTR","LTR/Copia","LTR/Gypsy","MobileElement","nonLTR","Other","Other/Simple","RC/Helitron","Retroelement","rRNA","Satellite","simple","SINE")
# add to reps
reps$rep_family <- repfam
reps$rep_family2 <- rf_c[repfam]
reps$Chr <- scaff2chrom[reps$scaff]
# get repeat coverage per-family in windows
reps.gr <- makeGRangesFromDataFrame(
  data.frame(chr = reps$Chr,
             start = reps$start,
             end = reps$end, 
             repeat_family = reps$rep_family2),
  keep.extra.columns = T
)
seqlevels(reps.gr) <- chr.info$chr.name
seqlengths(reps.gr) <- chr.info$length
reps.gr <- trim(reps.gr)
# make 1MB windows
bins <- GenomicRanges::tileGenome(seqlengths(reps.gr),
                                  tilewidth=1000000,
                                  cut.last.tile.in.chrom=TRUE)
# function to get repeat coverage
get_rep_coverage <- function(reps.gr, bins){
  fams <- sort(unique(reps.gr$repeat_family))
  rep_cov <- list()
  for(i in fams){
    # subset rep family
    my.reps <- reps.gr[which(reps.gr$repeat_family == i),]
    # get coverage
    my.cov <- coverage(my.reps)
    # set all regions covered by at least one rep to 1
    for(n in names(my.cov)){
      my.cov[[n]]@values <- ifelse(my.cov[[n]]@values > 0, 1, 0)
    }
    # get mean coverage across windows (i.e prop covered by a rep)
    my.rep_cov <- binnedAverage(bins = bins, numvar = my.cov, paste0(i, " coverage"))
    my.rep_cov <- as.data.frame(my.rep_cov)
    rep_cov[[i]] <- my.rep_cov
  }
  rep_cov <- Reduce(merge, rep_cov)
  colnames(rep_cov)[1] <- "Chr"
  rep_cov$start <- rep_cov$start-1 # so it matches bedtools windows
  rep_cov
}
rep_coverage <- get_rep_coverage(reps.gr = reps.gr, bins = bins)
# Marey map
mareymap <- Reduce(rbind, blocks_noOL)[,c("LG","scaffold", "Chr", "pos", "map.pos")]
# Recombination rate
rr <- Reduce(rbind, recomb)
rr$scaffold <- chrom2scaff[rr$Chr]
rr$start <- GC$start
rr$end <- GC$end
rr$mid <- rr$mid*1E6
# Combine windowed data
windows <- merge(NGenes, GC)
windows <- merge(windows, rr)
windows <- merge(windows, rep_coverage)
# plot correlations between stats
my.dat <- windows[,c("Chr", "N.genes", "pct_gc", "rr.lo", "Gypsy.coverage", "Copia.coverage")]
colnames(my.dat) <- c("chromosome","Genes/Mb", "GC %", "cM/Mb", "Gypsy density", "Copia density")
line_func <- function(data, mapping, ...) {
  ggplot(data, mapping) + 
    geom_point(size = 0.7) +
    geom_smooth(formula = y~x, method = lm, color = "black", se = TRUE,
                linetype = 1)
}
p <- ggpairs(my.dat, columns = 2:6, upper = list(continuous = wrap(ggally_cor, method = "spearman",size = 3.5)), lower=list(continuous=line_func), aes(colour=chromosome, alpha = 0.05), size = 0.8) + 
  theme_bw(base_size = 16) + 
  scale_color_manual(values = unname(chrom.cols)) 
p
ggsave("results/genome_struc_cor_all.pdf", width = 10, height = 10) 

#### 6. CIRCOS PLOT
# lengths
my.lengths <- chr.info$length
names(my.lengths) <- chr.info$chr.name
#genome
genome <- makeGRangesFromDataFrame(data.frame(chr = names(my.lengths), start=rep(1,12), end = my.lengths, strand=rep("*",12)))
seqlengths(genome) <- my.lengths
# marey map
mareymap.gr <- makeGRangesFromDataFrame(
  data.frame(chr = mareymap$Chr,
             start = mareymap$pos,
             end = mareymap$pos, 
             LG = mareymap$LG,
             map.pos = mareymap$map.pos),
  keep.extra.columns = T
)
seqlevels(mareymap.gr) <- names(my.lengths)
seqlengths(mareymap.gr) <- my.lengths
# recombination rate
windows.gr <- makeGRangesFromDataFrame(
  data.frame(chr = windows$Chr,
             start = windows$start,
             end = windows$end, 
             rr = (windows$rr.lo - min(windows$rr.lo, na.rm = T)) / sd(windows$rr.lo, na.rm = T),
             n.genes = (windows$N.genes - min(windows$N.genes)) / sd(windows$N.genes),
             GC = (windows$pct_gc - min(windows$pct_gc)) / sd(windows$pct_gc),
             Gypsy = (windows$Gypsy.coverage - min(windows$Gypsy.coverage)) / sd(windows$Gypsy.coverage),
             Copia = (windows$Copia.coverage - min(windows$Copia.coverage)) / sd(windows$Copia.coverage)),
  keep.extra.columns = T
)
seqlevels(windows.gr) <- names(my.lengths)
seqlengths(windows.gr) <- my.lengths
windows.gr <- trim(windows.gr)
# plot
ggbio() +
  circle(windows.gr, geom = 'bar', aes(y = GC), grid = F, fill = "darkgoldenrod",colour = "darkgoldenrod", size = 0.1) +
  circle(windows.gr, geom = 'bar', aes(y = Gypsy), grid = F, fill = "forestgreen",colour = "forestgreen", size = 0.1, ylim = c(0,1)) +
  circle(windows.gr, geom = 'bar', aes(y = Copia), grid = F, fill = "lightgreen",colour = "lightgreen", size = 0.1) +
  circle(windows.gr, geom = 'bar', aes(y = n.genes), size = 0.1, grid = F, fill = "steelblue",colour = "steelblue") +
  circle(subset(windows.gr, !is.na(rr)), geom = 'bar', aes(y = rr), size = 0.1, grid = F, fill = "lightcoral", colour = "lightcoral") +
  circle(mareymap.gr, geom = 'point', aes(y = map.pos), size = 0.1, grid = TRUE, grid.n = 4, grid.line = 'gray70', grid.background = "white", color = "black") + 
  circle(genome, geom = "ideo", aes(fill = seqnames)) + scale_fill_manual(values = chrom.cols) +
  circle(genome, geom = 'scale', size = 2) +
  circle(genome, geom = 'text', aes(label=seqnames), vjust=-1, size=3) +
  guides(fill="none")
ggsave("results/circos.pdf", width = 10, height = 10) 

#### 7. Repeat Manhattan plot
source("https://raw.githubusercontent.com/ogosborne/Scyliorhinus_canicula_popgen/refs/heads/main/code/manhattan2.R")
df <- windows
df$chrnum <- as.numeric(substr(df$Chr,4,5))
layout(matrix(c(rep(8,7),rep(c(1,2,3,4,5,6,7),4)), ncol = 5))
par(mar = c(1,6,0,1))
manhattan2(df,
           chr = "chrnum", 
           p = "Gypsy.coverage", 
           bp = "mid", 
           snp = "chrnum", 
           logp=FALSE, 
           ylab=expression(italic("Gypsy")), 
           cex.lab = 1.3, 
           chrlabs = rep("",12), 
           ylim = c(min(df$Gypsy.coverage), max(df$Gypsy.coverage)), 
           col = c("black", "cornflowerblue"), 
           las = 2,
           xlab = "",
           mgp=c(3.5,1,0), 
           xaxt = "n",
           type = "l",
           lwd = 2)
manhattan2(df,
           chr = "chrnum", 
           p = "Copia.coverage", 
           bp = "mid", 
           snp = "chrnum", 
           logp=FALSE, 
           ylab=expression(italic("Copia")), 
           cex.lab = 1.3, 
           chrlabs = rep("",12), 
           ylim = c(min(df$Copia.coverage), max(df$Copia.coverage)), 
           col = c("black", "cornflowerblue"), 
           las = 2,
           xlab = "",
           mgp=c(3.5,1,0), 
           xaxt = "n",
           type = "l",
           lwd = 2)
manhattan2(df,
           chr = "chrnum", 
           p = "Other.LTR.retrotransposon.coverage", 
           bp = "mid", 
           snp = "chrnum", 
           logp=FALSE, 
           ylab=expression(italic("Other LTR")), 
           cex.lab = 1.3, 
           chrlabs = rep("",12), 
           ylim = c(min(df$Other.LTR.retrotransposon.coverage), max(df$Other.LTR.retrotransposon.coverage)), 
           col = c("black", "cornflowerblue"), 
           las = 2,
           xlab = "",
           mgp=c(3.5,1,0), 
           xaxt = "n",
           type = "l",
           lwd = 2)
mtext("Proportion of nucleotides", side = 2, line = 7, cex = 2, at = 0)
manhattan2(df,
           chr = "chrnum", 
           p = "non.LTR.retrotransposon.coverage", 
           bp = "mid", 
           snp = "chrnum", 
           logp=FALSE, 
           ylab=expression(italic("Non-LTR\nRetrotransposon")), 
           cex.lab = 1.3, 
           chrlabs = rep("",12), 
           ylim = c(min(df$non.LTR.retrotransposon.coverage), max(df$non.LTR.retrotransposon.coverage)), 
           col = c("black", "cornflowerblue"), 
           las = 2,
           xlab = "",
           mgp=c(3.5,1,0), 
           xaxt = "n",
           type = "l",
           lwd = 2)
manhattan2(df,
           chr = "chrnum", 
           p = "DNA.transposon.coverage", 
           bp = "mid", 
           snp = "chrnum", 
           logp=FALSE, 
           ylab=expression(italic("DNA\nTransposon")), 
           cex.lab = 1.3, 
           chrlabs = rep("",12), 
           ylim = c(min(df$DNA.transposon.coverage), max(df$DNA.transposon.coverage)), 
           col = c("black", "cornflowerblue"), 
           las = 2,
           xlab = "",
           mgp=c(3.5,1,0), 
           xaxt = "n",
           type = "l",
           lwd = 2)
manhattan2(df,
           chr = "chrnum", 
           p = "Simple.repeat.coverage", 
           bp = "mid", 
           snp = "chrnum", 
           logp=FALSE, 
           ylab=expression(italic("Simple Repeat")), 
           cex.lab = 1.3, 
           chrlabs = chr.info$chr.name, 
           ylim = c( min(df$Simple.repeat.coverage), max(df$Simple.repeat.coverage)),
           col = c("black", "cornflowerblue"), 
           las = 1,
           mgp=c(3.5,1,0), 
           type = "l",
           lwd = 2)
