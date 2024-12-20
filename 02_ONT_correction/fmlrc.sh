# Illumina data
ill=/home/data/Illumina/F_read1and2_paired.fq.gz
# ONT data
ont=/home/data/ONT/MinION_all.fasta
# Correct
mkdir ./temp
gunzip -c ${ill} | awk 'NR % 4 == 2' | sort --parallel=40 -T ./temp | gzip > reads.sorted.txt.gz
rm -r ./temp
gunzip -c reads.sorted.txt.gz | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc-convert comp_msbwt.npy
fmlrc -p 40 comp_msbwt.npy ${ont} MinION_corr.fasta
