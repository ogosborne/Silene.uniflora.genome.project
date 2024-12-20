ont=MinION_corr.fasta
scaff=SLR_outdir/scaffold_set.fa

TGS-GapCloser.sh --scaff ${scaff} --reads ${ont} --output TGS_out --thread 20 --ne >TGS.log 2>TGS.err
 
