# train SNAP
fathom final_golden_genes.gff3.nr.golden.train.zff final_golden_genes.gff3.geneid.gff3.nr.golden.train.gff3.fasta -gene-stats | tee gene.statistics.log
fathom final_golden_genes.gff3.nr.golden.train.zff final_golden_genes.gff3.geneid.gff3.nr.golden.train.gff3.fasta -categorize 1000
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Pult . > Suniflora.hmm 
# set up predict dir
mkdir predict; cd predict
# split into scaffs
~/software/JAMg/bin/splitfasta.pl -i Su_softmasked_rn.fasta
# make cmd list
find Su_softmasked_rn.fasta_dir1 -maxdepth 1 -type f -exec sh -c \
 'echo "~/software/SNAP/snap ../train/Suniflora.hmm $1 -lcmask -quiet > $1.snap 2>/dev/null ; \
 ~/software/SNAP/zff2gff3.pl $1.snap > $1.snap.gff3 ; \
  ~/software/EVidenceModeler-1.1.1/EvmUtils/misc/SNAP_to_GFF3.pl $1.snap.gff3 > $1.snap.evm.gff3 "' \
  find-copy '{}' \; > snap.commands
# run SNAP prediction
~/software/JAMg/3rd_party/bin/ParaFly -shuffle -v -CPU 15 -c snap.commands -failed_cmds commands.list.failed | tee run.snap.log
# concatenate output
cat Su_softmasked.fasta_dir1/*.snap.evm.gff3 > snap.evm.gff3
# validate for EVM
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl snap.evm.gff3 

