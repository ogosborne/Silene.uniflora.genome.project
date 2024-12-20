ont=MinION_corr.fasta
scaff=Dovetail_v1.eq_mod.fa
# prepare the input bam files
bwa index ${scaff}
bwa mem -t 20 -a ${scaff} ${scaff} > align-self.sam
samtools view -Sb align-self.sam > align-self.bam
bwa mem -t 20 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y ${scaff} ${ont} > aligning.sam
samtools view -Sb aligning.sam > aligning.bam
# run SLR
SLR -c ${scaff} -r aligning.bam -d align-self.bam -p SLR_outdir -x 300 
