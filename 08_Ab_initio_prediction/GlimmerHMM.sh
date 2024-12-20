#### GlimmerHMM
# train: values from pasa post scripts and R (line 63)
trainGlimmerHMM \
 final_golden_genes.gff3.nr.golden.train.good.gb.fasta \
 final_golden_genes.gff3.nr.golden.train.good.gb.glimmer \
 -d Suniflora \
 -n 65425 \
 -f 143 \
 -l 309 >/dev/null
# set up predict dir
mkdir predict; cd predict
# split into scaffs
~/software/JAMg/bin/splitfasta.pl -i Su_hardmasked.fasta
# make cmd list
find Su_hardmasked.fasta_dir1 -maxdepth 1 -type f -exec sh -c \
 'echo "~/software/GlimmerHMM/bin/glimmerhmm_linux_x86_64 $1 ../train/Suniflora -f -g > $1.glimmer.gff 2>/dev/null ; \
  ~/software/EVidenceModeler-1.1.1/EvmUtils/misc/glimmerHMM_to_GFF3.pl $1.glimmer.gff > $1.glimmer.gff3 2>/dev/null ; \
  sed \"s/s_//g\" $1.glimmer.gff3 > $1.glimmer.gff3.rn 2>/dev/null ; \
  mv $1.glimmer.gff3.rn $1.glimmer.gff3"' \
  find-copy '{}' \; > glimmer.commands
# run 
~/software/JAMg/3rd_party/bin/ParaFly -shuffle -v -CPU 30 -c glimmer.commands -failed_cmds commands.list.failed | tee run.glimmer.log
# concatenate output
cat Su_hardmasked_.fasta_dir1/*.glimmer.gff3 > glimmer.evm.gff3
# validate
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl glimmer.evm.gff3 
