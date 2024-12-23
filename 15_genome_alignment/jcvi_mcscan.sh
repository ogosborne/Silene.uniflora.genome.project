### format data
# extract cds
for i in Suni Slat Scon ; do
   gffread -x ${i}.cds.fa -g ${i}.fasta ${i}.gff
done
# format fasta
for i in Suni Slat Scon ; do
   python -m jcvi.formats.fasta format ${i}.cds.fa ${i}.cds 
done
# format bed
for i in Suni Slat Scon ; do
   python -m jcvi.formats.gff bed --type=mRNA --key=Name ${i}.gff -o ${i}.bed
done

### MCScan
# find synteny
python -m jcvi.compara.catalog ortholog Slat Suni --cscore=.99 --no_strip_names --cpus=30 > Slat_Suni.out 2> Slat_Suni.err
python -m jcvi.compara.catalog ortholog Suni Scon --cscore=.99 --no_strip_names --cpus=30 > Suni_Scon.out 2> Suni_Scon.err
# make simple files
for i in Slat.Suni Suni.Scon ; do
   python -m jcvi.compara.synteny screen --minspan=30 --simple ${i}.anchors ${i}.anchors.simple
done
# make seqids file to order the chromosomes
echo -n '''Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,ChrX
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10''' > seqids
# add chromosome colours to anchor files
awk '{if ($3 ~ /Chr01/ ){print "firebrick*"$0 } 
else if ($3 ~ /Chr02/ ){print "red*"$0 } 
else if ($3 ~ /Chr03/ ){print "salmon*"$0 } 
else if ($3 ~ /Chr04/ ){print "darkorange*"$0 } 
else if ($3 ~ /Chr05/ ){print "gold*"$0 } 
else if ($3 ~ /Chr06/ ){print "yellowgreen*"$0 } 
else if ($3 ~ /Chr07/ ){print "forestgreen*"$0 } 
else if ($3 ~ /Chr08/ ){print "teal*"$0 } 
else if ($3 ~ /Chr09/ ){print "steelblue*"$0 } 
else if ($3 ~ /Chr10/ ){print "cornflowerblue*"$0 } 
else if ($3 ~ /Chr11/ ){print "slateblue*"$0 } 
else if ($3 ~ /Chr12/ ){print "indigo*"$0 } }' Slat.Suni.anchors.simple > Slat.Suni.anchors.simple.cols
awk '{if ($1 ~ /Chr01/ ){print "firebrick*"$0 } 
else if ($1 ~ /Chr02/ ){print "red*"$0 } 
else if ($1 ~ /Chr03/ ){print "salmon*"$0 } 
else if ($1 ~ /Chr04/ ){print "darkorange*"$0 } 
else if ($1 ~ /Chr05/ ){print "gold*"$0 } 
else if ($1 ~ /Chr06/ ){print "yellowgreen*"$0 } 
else if ($1 ~ /Chr07/ ){print "forestgreen*"$0 } 
else if ($1 ~ /Chr08/ ){print "teal*"$0 } 
else if ($1 ~ /Chr09/ ){print "steelblue*"$0 } 
else if ($1 ~ /Chr10/ ){print "cornflowerblue*"$0 } 
else if ($1 ~ /Chr11/ ){print "slateblue*"$0 } 
else if ($1 ~ /Chr12/ ){print "indigo*"$0 } }' Suni.Scon.anchors.simple > Suni.Scon.anchors.simple.cols
# make layout file
echo -n '''# y, xstart, xend, rotation, color, label, va,  bed
 .9,     .25,    .8,       0,     k, S. lat, top, Slat.bed
 .7,     .25,    .8,       0,     k, S. uni, top, Suni.bed
 .5,     .25,    .8,       0,     k, S. con, top, Scon.bed
# edges
e, 0, 1, Slat.Suni.anchors.simple.cols
e, 1, 2, Suni.Scon.anchors.simple.cols''' > layout
# make plot
python -m jcvi.graphics.karyotype --font=Arial -o karyotype.pdf seqids layout 
