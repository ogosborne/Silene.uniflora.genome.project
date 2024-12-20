
# run busco for Augustus training parameters
genome=Su_softmasked.fasta
export BUSCO_CONFIG_FILE=/home/software/busco/config/myconfig.ini
export AUGUSTUS_CONFIG_PATH="/home/software/Augustus/config/"
export PATH="/home/software/Augustus/bin:$PATH"
export PATH="/home/software/Augustus/scripts:$PATH"
busco -m genome --augustus_species tomato -i ${genome} -o Su -l eudicots_odb10 --cpu 20
# get hints
cat final_golden_genes.gff3.nr.golden.optimization.good.gb.fasta final_golden_genes.gff3.nr.golden.train.good.gb.fasta > complete_reference_training_set.fasta
# prepare hints for training set reference: RNA-Seq and repeatmasking/exons
~/software/JAMg/bin/align_rnaseq_gsnap.pl -fasta complete_reference_training_set.fasta -gmap_dir . -suffix -dbname reference_training -pattern1 _1 -pattern2 _2 -nofail
~/software/JAMg/bin/augustus_RNAseq_hints.pl -bam *uniq_mult.bam -genome complete_reference_training_set.fasta

# run
augustus --species=Su --singlestrand=false  --extrinsicCfgFile=extrinsic.E.cfg --noInFrameStop=true --UTR=off --hintsfile=uniq_mult.bam.hints $genome > augustus.busco.gff
