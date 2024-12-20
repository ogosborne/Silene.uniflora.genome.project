######################### PREPARE INPUTS

#### Stringtie
# validate gff3
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl stringtie.transdecoder.genome.gff3

#### PASA 
# rename transdecoder to pasa in pasa gff
sed "s/transdecoder/pasa/g" mydb_pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 > pasa.transdecoder.genome.gff3
# validate 
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl pasa.transdecoder.genome.gff3

#### Exonerate
# convert
for i in *.exonerate.gff ; do
	Exonerate_to_evm_gff3.pl $i > ${i}.gff3
done
cat *.exonerate.gff.gff3 > Exonerate.all.gff3
# validate
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl Exonerate.all.gff3

#### GeneMark-ES
# convert
~/software/EVidenceModeler-1.1.1/EvmUtils/misc/GeneMarkHMM_GTF_to_EVM_GFF3.pl genemark.gtf > genemark.gff3
# validate 
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl genemark.gff3

#### Augustus
# convert
~/software/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl augustus.busco.gff > augustus.gff3
# validate
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl augustus.gff3

#### GlimmerHMM 
# convert
~/software/EVidenceModeler-1.1.1/EvmUtils/misc/glimmerHMM_to_GFF3.pl Glimmer_output.gff > glimmer.evm.gff3
# validate gff3
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl Glimmer_output.gff3

#### SNAP
# convert
~/software/EVidenceModeler-1.1.1/EvmUtils/misc/SNAP_to_GFF3.pl SNAP_output.gff3 > snap.evm.gff3
# validate
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl snap.evm.gff3

# combine ab initio and RNA gene predictions
cat augustus.gff3 genemark.gff3 glimmer.evm.gff3 snap.evm.gff3 stringtie.transdecoder.genome.gff3 pasa.transdecoder.genome.gff3 > gene.predictions.gff3 

#### make weights file

# Exonerate
echo -e "PROTEIN\texonerate\t5" > weights.txt
# PASA 
echo -e "OTHER_PREDICTION\tpasa\t10" >> weights.txt
# stringtie
echo -e "OTHER_PREDICTION\tstringtie\t10" >> weights.txt
# Augustus: higher than other ab initio as it includes extrinsic evidence
echo -e "ABINITIO_PREDICTION\tAugustus\t2" >> weights.txt
# Genemark
echo -e "ABINITIO_PREDICTION\tGeneMark.hmm\t1" >> weights.txt
# SNAP 
echo -e "ABINITIO_PREDICTION\tSNAP\t1" >> weights.txt
# GlimmerHMM
echo -e "ABINITIO_PREDICTION\tGlimmerHMM\t1" >> weights.txt

######################### RUN

#### partition inputs
~/software/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl \
--genome Su_softmasked.fasta \
--gene_predictions gene.predictions.gff3 \
--protein_alignments Exonerate.all.gff3 \
--segmentSize 1000000 \
--overlapSize 100000 \
--partition_listing partitions.txt

#### write command list
~/software/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl \
--genome Su_softmasked.fasta \
--weights weights.txt \
--gene_predictions gene.predictions.gff3 \
--protein_alignments Exonerate.all.gff3 \
--output_file_name evm.out \
--partitions partitions.txt > commands.list

#### run EVM using parafly
~/software/JAMg/3rd_party/bin/ParaFly -shuffle -v -CPU 35 -c commands.list -failed_cmds commands.list.failed

#### recombine partitions
~/software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions.txt --output_file_name evm.out

#### convert to gff3
# convert
~/software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions.txt --output evm.out --genome Su_softmasked.fasta
# combine 
find . -name evm.out.gff3 -exec cat '{}' \; >> all.evm.gff3
