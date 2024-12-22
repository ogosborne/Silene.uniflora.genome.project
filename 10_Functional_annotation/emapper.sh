
gffread -g Su_softmasked.fasta -S -y PASA.final.gff3.pep.fasta PASA.final.gff3

emapper.py \
--tax_scope 33090 \
--output_dir results \
--data_dir data \
--output eggNOG_mapper \
--query_cover 70 \
--subject_cover 70 \
--evalue 0.0001 \
-m diamond \
--pfam_realign realign \
--cpu 20 \
--num_servers 20 \
-i PASA.final.gff3.pep.fasta | tee emapper.log

