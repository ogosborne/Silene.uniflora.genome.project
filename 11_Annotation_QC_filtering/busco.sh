for lineage in viridiplantae_odb10 embryophyta_odb10 eudicots_odb10 ; do
	busco -m protein -i PASA.final.gff3.pep.fasta -o PASA.final_${lineage} -l ${lineage} --cpu 20 -f --offline
	cp ./PASA.final_${lineage}/run_${lineage}/full_table.tsv PASA.final_${lineage}_full_table.tsv
done
