#!/usr/bin/env bash

conform="/mnt/lustre/groups/CBBI1243/KEVIN/bioTools/conform-gt.24May16.cee.jar"
plink2="${HOME}/bin/plink2"
#-- Convert PLINK binary to VCF (only common SNPs MAF>0.01)
#for chr in {1..22} X; do
#$plink2 \
#        --bfile /mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/analysis/qc-camgwas \
#        --chr $chr \
#	--maf 0.01 \
#	--geno 0.05 \
#        --export vcf-4.2 id-paste=fid id-delim="_" bgz \
#        --real-ref-alleles \
#        --out /mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/analysis/chr${chr}
#
#        #java -jar ${conform} ref=/mnt/lustre/groups/1000genomes/annotation/REF_IMPUTATION/refdata/chr${chr}.1kg.phase3.v5a.vcf.gz gt=chr${chr}.vcf.gz chrom=${chr} out=new.chr${chr} excludesamples=conform.list.excl match=POS
#
#done

for i in {1..22}; do echo $i; done | parallel echo "-f -p vcf /mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/analysis/chr{}.vcf.gz" | xargs -P10 -n4 tabix

#-- Check SNP against 1KGP3 for overlap and strand using BEAGLE conform script
for i in {1..22}; do echo $i; done | parallel echo "-jar ${conform} ref=/mnt/lustre/groups/1000genomes/annotation/REF_IMPUTATION/refdata/chr{}.1kg.phase3.v5a.vcf.gz gt=/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/analysis/chr{}.vcf.gz chrom={} out=/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/new/new.chr{} match=POS excludesamples=/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/analysis/conform.list.excl" | xargs -P5 -n8 java

#-- Merge the new files with allele overlap in ref, removing dups then separate again removing multi-allelic
for i in {1..22}; do echo $i; done | parallel echo "-f -p vcf new.chr{}.vcf.gz" | xargs -P10 -n4 tabix
bcftools concat -a -D --threads 24 -Oz new.chr{1..22}.vcf.gz -o camgwas.vcf.gz
bcftools index --threads 24 -ft camgwas.vcf.gz
for i in {1..22}; do echo $i; done | parallel echo "view -m2 -M2 -r {} -Oz -o chr{}.vcf.gz camgwas.vcf.gz" | xargs -P10 -n9 bcftools

#-- Phase with BEAGLEv5.1 without imputation
#for i in {1..22} X; do echo $i; done | parallel echo "gt=/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/new/chr{}.vcf.gz ref=/mnt/lustre/groups/CBBI1243/KEVIN/db/Phase3_merged.vcf.gz chrom={} nthreads=5 impute=false ne=20000 out=/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/new/chr{}-phased" | xargs -P5 -n7 $beagle5

#-- Phase while imputing
for i in {1..22}; do echo $i; done | parallel echo "gt=/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/new/chr{}.vcf.gz ref=/mnt/lustre/groups/1000genomes/annotation/REF_IMPUTATION/refdata/chr{}.1kg.phase3.v5a.vcf.gz chrom={} impute=true ne=20000 gp=true out=/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/new/chr{}-imputed" | xargs -P5 -n7 $beagle5

#-- Merge imputed
for i in {1..22}; do echo $i; done | parallel echo "-f -p vcf chr{}-imputed.vcf.gz" | xargs -P10 -n4 tabix
bcftools concat -a -D --threads 24 -Oz -o imputed.vcf.gz chr{1..22}-imputed.vcf.gz
bcftools index --threads 24 -ft imputed.vcf.gz
bcftools view -v snps,indels -m2 -M2 -i 'DR2>0.75' imputed.vcf.gz --threads 24 -Oz -o filtered_r20.75_imputed.vcf.gz
bcftools index --threads 24 -ft filtered_r20.75_imputed.vcf.gz
