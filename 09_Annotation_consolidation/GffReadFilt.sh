#### get lists of genes which pass gffread filters
# multi-exon mRNAs that have any intron with a non-canonical splice site consensus
~/software/gffread/gffread --keep-genes -N --stream -g Su_softmasked.fasta --keep-genes PASA.final.gff3 | awk '{if ($3 == "mRNA") print $0;}' | cut -f9 | cut -d";" -f2 |  sed 's/geneID=//g' > IDs.pass.N.txt
# mRNAs with CDS having in-frame stop codons 
~/software/gffread/gffread --keep-genes -V --stream -g Su_softmasked.fasta --keep-genes PASA.final.gff3 | awk '{if ($3 == "mRNA") print $0;}' | cut -f9 | cut -d";" -f2 |  sed 's/geneID=//g' > IDs.pass.V.txt
# single exon transcripts
~/software/gffread/gffread --keep-genes -U --stream -g Su_softmasked.fasta --keep-genes PASA.final.gff3 | awk '{if ($3 == "mRNA") print $0;}' | cut -f9 | cut -d";" -f2 |  sed 's/geneID=//g' > IDs.pass.U.txt
# mRNAs that either lack initial START codon or the terminal STOP codon, or have an in-frame stop codon
~/software/gffread/gffread --keep-genes -J --stream -g Su_softmasked.fasta --keep-genes PASA.final.gff3 | awk '{if ($3 == "mRNA") print $0;}' | cut -f9 | cut -d";" -f2 |  sed 's/geneID=//g' > IDs.pass.J.txt
