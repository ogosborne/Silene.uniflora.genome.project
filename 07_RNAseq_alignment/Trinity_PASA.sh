# Su_softmasked.fasta is the final assembly following repeat masking by RepeatMasker (all repeats converted to lowercase)

# run Trinity
Trinity --seqType fq --left Cleaned_Root_Flower_Leaf_1.fq --right Cleaned_Root_Flower_Leaf_2.fq --CPU 40 --max_memory 300G 

### run seqclean
$PASAHOME/bin/seqclean Trinity.fasta

### run accession_extractor.pl    
$PASAHOME/misc_utilities/accession_extractor.pl < Trinity.fasta > tdn.accs

### align transcripts to genome with PASA
$PASAHOME/Launch_PASA_pipeline.pl -c sqlite.alignAssembly.config -C -R -g Su_softmasked.fasta --ALIGNERS blat,gmap --TRANSDECODER --CPU 20 -T -t Trinity.fasta.clean -u Trinity.fasta --TDN tdn.accs |& tee pasa.log

### run transdecoder on PASA to get genome-coordinate gff3
$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta mydb_pasa.sqlite.assemblies.fasta --pasa_transcripts_gff3 mydb_pasa.sqlite.pasa_assemblies.gff3 |& tee pasa_asmbls_to_training_set.log
