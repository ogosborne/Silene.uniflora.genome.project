# trim data
for s in Root Flower Leaf ; do
	java -jar trimmomatic-0.39.jar PE -threads 3 -phred33 \
	-trimlog ${s}_trimlog.txt \
	-summary ${s}_summary.txt \
	${s}_R1_001.fq.gz \
	${s}_R2_001.fq.gz \
	${s}_trimmed_R1.fq.gz \
	${s}_trimmed_U1.fq.gz \
	${s}_trimmed_R2.fq.gz \
	${s}_trimmed_U2.fq.gz \
	ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:TRUE \
	LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:70
done
# concatenate data
cat *_trimmed_R1.fq.gz | zcat > Cleaned_Root_Flower_Leaf_1.fq
cat *_trimmed_R2.fq.gz | zcat > Cleaned_Root_Flower_Leaf_2.fq
