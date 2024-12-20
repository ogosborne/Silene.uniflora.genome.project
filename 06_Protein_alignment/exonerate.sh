# Su_hardmasked.fasta is the final assembly following repeat masking by RepeatMasker (all repeats converted to N)
# Amaranthus hypochondriacus
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Ahypochondriacus_459_v2.1.protein.fa --target Su_hardmasked.fasta > Ah.out
# Arabidopsis thaliana
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Athaliana_447_Araport11.protein.fa --target Su_hardmasked.fasta > At.out
# Daucus carota
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Dcarota_388_v2.0.protein.fa  --target Su_hardmasked.fasta > Dc.out
# Helianthus annuus
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Hannuus_494_r1.2.protein.fa  --target Su_hardmasked.fasta > Ha.out
# Lactuca sativa
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Lsativa_467_v5.protein.fa  --target Su_hardmasked.fasta > Lc.out
# Mimulus guttatus
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Mguttatus_256_v2.0.protein.fa  --target Su_hardmasked.fasta > Mg.out
# Olea europaea
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Oeuropaea_451_v1.0.protein.fa  --target Su_hardmasked.fasta > Oe.out
# Solanum lycopersicum
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Slycopersicum_514_ITAG3.2.protein.fa --target Su_hardmasked.fasta > Sl.out
# Solanum tuberosum
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" --query Stuberosum_448_v4.03.protein.fa --target Su_hardmasked.fasta > St.out
