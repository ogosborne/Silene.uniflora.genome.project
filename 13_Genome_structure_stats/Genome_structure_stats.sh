# bedtools version: 2.27.1
# bedops version (gff2bed): 2.4.41
# dirs
mkdir results/gc
mkdir results/gene_density
mkdir results/windows
# env
conda activate bedops
fai=Su_softmasked.fasta.fai
gen=Su_softmasked.fasta
gff=annotation/best.filt.annots.gff3
# Make bed file with 1Mb windows
sort -k2,2n $fai | tail -n12 | cut -f1,2 | tac > results/windows/genome.tsv
bedtools makewindows -g results/windows/genome.tsv -w 1000000 > results/windows/1Mb_window.bed
# gc content 
bedtools nuc -fi $gen -bed results/windows/1Mb_window.bed > results/gc/GC_1Mb.tsv
# gene density
awk '$3=="gene" {print $0}' $gff | gff2bed > results/gene_density/genes.bed
bedtools intersect -a results/windows/1Mb_window.bed -b results/gene_density/genes.bed -c > results/gene_density/ngenes_1Mb.tsv 
