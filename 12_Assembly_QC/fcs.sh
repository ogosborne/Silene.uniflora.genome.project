# screen genome for adaptors
genome=Su_softmasked.fasta 
run_fcsadaptor.sh --fasta-input $genome --output-dir ./results --euk --container-engine singularity --image fcs-adaptor.sif
# modify output to fix (mask) adaptors rather than splitting contigs
sed 's/ACTION_TRIM/FIX/g' results/fcs_adaptor_report.txt > results/fcs_adaptor_report_fix.txt
# mask adaptors
cat $genome | python3 software/fcs.py clean genome --action-report results/fcs_adaptor_report_fix.txt --output results/Su_softmasked_maskAdaptors.fasta --contam-fasta-out results/contam.fasta
# screen genome for contaminant organisms
taxid=39919
genome=results/Su_softmasked_maskAdaptors.fasta
python3 software/fcs.py screen genome --fasta $genome2 --out-dir results2 --gx-db $LOCAL_DB --tax-id $taxid 

