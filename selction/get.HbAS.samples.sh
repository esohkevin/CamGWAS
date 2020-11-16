#!/usr/bin/env bash

base="/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/"
new="${base}new/"
sel="${new}selction/"

bcftools view -i 'GT=="het"' -r 11:5248232 ${new}filtered_r20.75_imputed.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' | cut -f5- | sed 's/\t/\n/g' | sed 's/=/\t/g' | awk '$2=="0|1" || $2=="1|0" {print $1"\t"$1"\t"$2}' > ${sel}HbAS.samples.txt
bcftools view -i 'GT=="het"' -r 11:5248232 ${new}filtered_r20.75_imputed.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' | cut -f5- | sed 's/\t/\n/g' | sed 's/=/\t/g' | awk '$2=="0|0" {print $1"\t"$1"\t"$2}' > ${sel}HbAA.samples.txt
bcftools view -i 'GT=="het"' -r 11:5248232 ${new}filtered_r20.75_imputed.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' | cut -f5- | sed 's/\t/\n/g' | sed 's/=/\t/g' | awk '$2=="1|1" {print $1"\t"$1"\t"$2}' > ${sel}HbSS.samples.txt

cat HbAS.samples.txt HbSS.samples.txt | cut -f1 > haps.HbS.chroms.samples.txt

bcftools view -v snps -m2 -M2 -i 'AF>=0.01' -r 11:4248232-6248232 ${new}filtered_r20.75_imputed.vcf.gz --threads 24 -Oz -o ${sel}hbb.vcf.gz

bcftools view -S haps.HbS.chroms.samples.txt --force-samples -r  11:5291563,11:5269799,11:5263683,11:5260458 ${new}filtered_r20.75_imputed.vcf.gz --threads 24 -Oz -o ${sel}haps.vcf.gz

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GT]\n' haps.vcf.gz > HbS.chroms.samples.haps.txt

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' haps.vcf.gz | sed 's/|/,/g' | sed 's/\t/,/g' > haps.csv

bcftools query -S ba.haps.HbS.chroms.samples.txt -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' haps.vcf.gz | sed 's/|/,/g' | sed 's/\t/,/g' > ba.haps.csv
bcftools query -S sb.haps.HbS.chroms.samples.txt -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' haps.vcf.gz | sed 's/|/,/g' | sed 's/\t/,/g' > sb.haps.csv

#   Haplotype rs3834466 rs28440105 rs10128556 rs968857 
#   AI        GT        C          T          T
#   SEN       G         C          T          T
#   BEN       G         C          C          T
#   CAR       G         C          C          C
#   CAM       G         A          C          T

#   HBB HinfI: 11:5246356 rs10837631
#   HBB AvaII: 11:5247791 rs10768683
#   HBP1 HincII: rs968857
#   HBP1 HincII: 11:5263683 rs10128556
#   HBG1 HincII: rs28440105
#   HBG2 HindIII: rs2070972
#   HBG2 XmnI: rs7482144
#   HBE1 HincII: rs3834466 


#   HBBP1 HincII: 11:5260458 rs968857
#   HBG1 (gamma A) HindIII: 11:5269799 rs28440105
#   HBG2 (gamma G) HindIII: 11:5274717 rs2070972
#   HBE1 HincII: 11:5291563 rs3834466 
