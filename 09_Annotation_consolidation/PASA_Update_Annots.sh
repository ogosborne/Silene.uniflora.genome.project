
### after EVM, rerun PASA to add UTR and alternatively spliced transcripts

# Load EVM gene annotations
$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi \
     -c sqlite.alignAssembly.config \
     -g Su_softmasked.fasta \
     -P all.evm.gff3 |& tee Load_Current_Gene_Annotations1.log

# Update EVM annotations
$PASAHOME/Launch_PASA_pipeline.pl \
        -c annotCompare.config \
        -A \
        -g Su_softmasked.fasta \
        -t Trinity.fasta.clean |& tee pasa.annotCompare1.log
# Trinity.fasta.clean created in 07_RNAseq_alignmen/Trinity_PASA.sh

# load new PASA annotations (1)
$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi \
     -c sqlite.alignAssembly.config \
     -g Su_softmasked.fasta \
     -P mydb_pasa.sqlite.gene_structures_post_PASA_updates.27063.gff3 |& tee Load_Current_Gene_Annotations2.log

# Update PASA annotations (1)
$PASAHOME/Launch_PASA_pipeline.pl \
        -c annotCompare.config \
        -A \
        -g Su_softmasked.fasta \
        -t Trinity.fasta.clean |& tee pasa.annotCompare2.log
mv 4.mydb_pasa.sqlite.gene_structures_post_PASA_updates.25334 PASA.final.gff3

# load new PASA annotations (2)
$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi \
     -c /sqlite.alignAssembly.config \
     -g Su_softmasked.fasta \
     -P mydb_pasa.sqlite.gene_structures_post_PASA_updates.34512.gff3 |& tee Load_Current_Gene_Annotations3.log

# Update PASA annotations (2)
$PASAHOME/Launch_PASA_pipeline.pl \
        -c annotCompare.config \
        -A \
        -g Su_softmasked.fasta \
        -t Trinity.fasta.clean |& tee pasa.annotCompare3.log
