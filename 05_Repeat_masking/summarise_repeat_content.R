library(dplyr)
# Repeat content
reps <- read.table("data/RepeatMasker/Su_DT_SLR_TGSGC_repeats.gff", sep = "\t")[,c(1,4,5,9)]
colnames(reps) <- c("scaff", "start", "end", "info")
reps$scaff <- paste0("s_", reps$scaff)
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
# length
reps$length <- reps$end - reps$start + 1
# Total repeat bases
total_repeat_bases <- sum(reps$length)
# Total genome size 
genome_size <- 1269363085  
# Overall repeat proportion
repeat_percent <- (total_repeat_bases / genome_size) * 100
repeat_percent
# Repeat proportion by family
repfamily_summary <- reps %>%
  mutate(length = end - start + 1) %>%
  group_by(rep_family2) %>%
  summarise(
    n_elements = n(),
    total_bases = sum(length),
    percent_of_genome = (sum(length) / genome_size) * 100
  ) %>%
  arrange(desc(percent_of_genome))

repfamily_summary
