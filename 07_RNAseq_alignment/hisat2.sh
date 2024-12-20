genome=Su_hardmasked.fasta
r1=Cleaned_Root_Flower_Leaf_1.fq
r2=Cleaned_Root_Flower_Leaf_2.fq
# build index
hisat2-build ${genome} genome
# map
hisat2 -x genome -1 ${r1} -2 ${r2} -S mapping.sam -p 40 --dta  2>&1 | tee hisat2.stat 
# sort
samtools sort -T . -m 8G -@ 20  mapping.sam > sorted.bam
