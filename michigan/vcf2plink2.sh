#!/usr/bin/env bash


plink2="${HOME}/bin/plink2"
home="/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/"
phase="${home}phase/"
sb="${phase}sbanimp/"
ba="${phase}banimp"
fm="${phase}5M/"
an="${home}analysis/"
as="${home}assoc_results/"
asd="${an}assoc/"
hm="/home/kesoh/bin/"
em="${an}emmax/"
gc="${an}gcta/"
gcta="${HOME}/bin/gcta64"
mc="/mnt/lustre/groups/CBBI1243/KEVIN/git/GWAS/new/mich/"

#plink2 --vcf dosage/chr22.dose.vcf.gz 'dosage=DS' --import-dosage-certainty 0.90 --extract-col-cond dosage/chr22.info.gz 7 1 0 --vcf-require-gt --maf 0.01 --geno 0.05 --out chr22 --make-bed --extract-col-cond-min 0.75

#for i in {1..22}; do echo $i; done | parallel echo "--vcf dosage/chr{}.dose.vcf.gz 'dosage=GP' --import-dosage-certainty 0.90 --pheno ${an}raw-camgwas.fam --pheno-col-nums 6 --update-sex ${an}raw-camgwas.fam col-num=5 --make-bed --double-id --out chr{} --extract-col-cond dosage/chr{}.info.gz 7 1 0 --extract-col-cond-min 0.65 --max-alleles 2 --snps-only just-acgt" | xargs -I input -P5 sh -c "plink2 input"
#for i in {1..22}; do echo chr$i; done > merge.list
#plink --merge-list merge.list --keep-allele-order --out camgwas
  
# plink2 --bfile camgwas --maf 0.01 --mind 0.10 --geno 0.10 --hwe 1e-08 --make-bed --out clean --keep ${asd}sbba.pca.glm.cov

#---------------------------------------------------------
#for i in {1..22}; do echo $i; done | parallel echo "--bfile cam.clean --pca 50 --chr {} --out chr{}.PC50" | xargs -I input -P5 sh -c "plink2 input"

for p in sb; do for i in {1..22}; do echo $i; done | parallel echo "--bfile clean --glm sex hide-covar --condition 11:5248232:T:A --adjust --ci 0.95 --out ${p}.chr{}.assoc --chr {} --covar ${asd}sbba.pca.glm.cov --covar-name PC1-PC20 --keep ${asd}${p}.pca.glm.cov" | xargs -I input -P5 sh -c "plink2 input"; done

for m in hethom dominant recessive genotypic; do for p in sbba; do for i in {1..22}; do echo $i; done | parallel echo "--bfile clean --glm sex hide-covar ${m} --condition 11:5248232:T:A --adjust --ci 0.95 --out ${p}.${m}.chr{}.assoc --chr {} --covar ${asd}sbba.pca.glm.cov --covar-name INTERCEPT,PC1-PC20 --keep ${asd}${p}.pca.glm.cov" | xargs -I input -P5 sh -c "plink2 input"; done; done

for k in sb; do for j in dominant recessive genotypic hethom; do cat ${k}.${j}.chr{1..22}.assoc.PHENO1.glm.logistic | sed '1d' | awk '{print $1,$2,$3,$4"/"$5,$7,$8,$9,$11,$12,$14}' | sed '1 i CHR POS SNP REF/ALT TEST OBS OR L95 U95 P' | grep -v -e NA -e CHROM | sort -g -k10 > ${k}.${j}.assoc; done; done

for k in sb; do cat ${k}.chr{1..22}.assoc.PHENO1.glm.logistic | sed '1d' | awk '{print $1,$2,$3,$4"/"$5,$7,$8,$9,$11,$12,$14}' | sed '1 i CHR POS SNP REF/ALT TEST OBS OR L95 U95 P' | grep -v -e NA -e CHROM | sort -g -k10 > ${k}.assoc; done

for i in *.assoc; do echo $i; done | parallel echo ${an}fdr.R {} | xargs -I input -P10 sh -c "Rscript input"
for i in *.assoc.adj.txt; do mv $i ${i/.adj.txt/}; done
#for i in *.plot.gz; do echo $i; done | parallel echo ${an}assoc.R {} | xargs -I input -P 10 sh -c "Rscript input"

#for i in sb; do for j in hethom recessive dominant; do cat ${i}.${j}.chr{1..22}.assoc.PHENO1.glm.logistic.adjusted | cut -f2,4,9 | grep -v "UNADJ" | sed 's/:/\t/g' | cut -f1-2,5-6 | sort -g -k3 | sed 's/\t/ /g' | sed '1 i CHR POS P P_BH' | gzip -c > ${i}.${j}.plot.gz; done; done

#------------------------------------------------------
# for i in {1..22}; do echo $i; done | parallel echo "--bfile cam.clean --pca 50 --chr {} --out chr{}.PC50" | xargs -I input -P5 sh -c "plink2 input"
# cp ${asd}*.cov .
# for k in *.cov; do for j in dominant recessive genotypic hethom; do for i in {1..22}; do echo $i; done | parallel echo "--bfile cam.clean --glm sex hide-covar ${j} --adjust --out chr{}.${j}.${k/.cov/} --chr {} --ci 0.95 --covar chr{}.PC50.eigenvec --covar-name PC1-PC50 --keep ${k}" | xargs -I input -P5 sh -c "plink2 input"; done; done
# 
# for k in *.cov; do for j in dominant recessive genotypic hethom; do cat chr{1..22}.${j}.${k/.cov/}.PHENO1.glm.logistic | sed '1d' | awk '{print $1,$2,$3,$4"/"$5,$7,$8,$9,$11,$12,$14}' | sed '1 i CHR POS SNP REF/ALT TEST OBS OR L95 U95 P' | grep -v -e NA -e CHROM | sort -g -k10 > ${j}.${k/.cov/}.assoc; done; done

#for i in *.assoc; do echo $i; done | parallel echo ${an}fdr.R {} | xargs -I input -P5 sh -c "Rscript input"
#for i in *.assoc.adj.txt; do mv $i ${i/.adj.txt/}; done
#for i in *.assoc; do echo $i; done | parallel echo ${an}assoc.R {} | xargs -I input -P5 sh -c "Rscript input"
