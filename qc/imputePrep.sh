#!/usr/bin/env bash

kgp="../phase/1000GP_Phase3/"
eig="${baseDir}/popstruct/eig/"
samples="../samples/"
k="../1000G/"
conform="/mnt/lustre/groups/CBBI1243/KEVIN/bioTools/conform-gt.24May16.cee.jar"
plink2="${HOME}/bin/plink2"
#echo -e """\e[1m\e[38;5;40m
# #	 Extracting eigen-corrected samples
# #	------------------------------------
# #	\e[0m
# #"""
# #plink \
# #        --bfile qc-camgwas \
# #        --keep ${samples}eig.ids \
# #        --make-bed \
# #	--autosome \
# #	--keep-allele-order \
# #        --out qc-camgwas-eig-cor
# #
# #echo -e """\e[1m\e[38;5;40m 
# #	 Getting Palimcromic SNPs
# #	--------------------------
# #	\e[0m
# #"""
# 
# # Compute allel frequencies
plink \
	--bfile qc-camgwas \
	--allow-no-sex \
	--freq \
	--out qc-camgwas

# echo -e """\e[1m\e[38;5;40m
# 	 Extracting palindromic SNPs in R for subsequent exclusion
# 	------------------------------------------------------------
# 	\e[0m
# """
Rscript palindromicsnps.R qc-camgwas.bim
# 
# echo -e """\e[1m\e[38;5;40m
# 	 Checking Strand using the Wrayner perl script
# 	------------------------------------------------
# 	\e[0m
# """
# ./checkstrand.sh
# 
# # Add command line to remove at-cg SNPs
# plink --bfile qc-camgwas --exclude at-cg.snps --make-bed --out TEMP0
# plink --bfile TEMP0 --allow-no-sex --exclude Exclude-qc-camgwas-1000G.txt --make-bed --out TEMP1
# plink --bfile TEMP1 --update-chr Chromosome-qc-camgwas-1000G.txt 2 1 --make-bed --out TEMP2
# plink --bfile TEMP2 --update-map Position-qc-camgwas-1000G.txt 2 1 --make-bed --out TEMP3
# plink --bfile TEMP3 --flip Strand-Flip-qc-camgwas-1000G.txt --make-bed --out qc-camgwas-updated
# #plink2 --bfile TEMP4 --ref-allele force refSites.txt 4 1 --make-bed --out qc-camgwas-updated
# #plink2 --bfile TEMP4 --ref-allele force chimp_anc.txt 3 2 --make-bed --out TEMP5
# #plink2 --bfile TEMP5 --ref-allele force ancSites.txt 4 3 --make-bed --out qc-camgwas-updated
# 
# echo -e """\e[1m\e[38;5;40m
# 	 Converting plink binary files to VCF
# 	--------------------------------------
# 	\e[0m
# """
# plink2 \
#         --bfile qc-camgwas-updated \
#         --chr $chr
#         --ref-allele force ancSites.txt 4 2 \
#         --export vcf-4.2 id-paste=fid bgz \
#         --real-ref-alleles \
#         --out qc-camgwas-updated

# tabix -f -p vcf qc-camgwas-updated.vcf.gz
# bcftools sort qc-camgwas-updated.vcf.gz -Oz -o qc-camgwas-updated.vcf.gz
# bcftools index --tbi qc-camgwas-updated.vcf.gz
# 
# #mv qc-camgwas-updated.vcf.gz ../phase/
# 
# #./plink2vcf.sh
# 
# #echo -e """
# #	 Chechk VCF for errors using the checkVCF.py script
# #	----------------------------------------------------
# #"""
# 
# #----Activate Python2.7 env for checkVCF
# source ~/.bashrc
# conda activate py2
# 
# ./vcfcheck.sh
# 
# conda deactivate
# 
# 
# rm TEMP*


#conform="/mnt/lustre/groups/CBBI1243/KEVIN/bioTools/conform-gt.24May16.cee.jar"

#for chr in {1..22} X; do
#    plink2 \
#        --bfile qc-camgwas \
#        --chr $chr \
#        --export vcf-4.2 id-paste=fid bgz \
#        --real-ref-alleles \
#        --out chr${chr}
#done

# for i in {1..22} X; do bcftools index -ft --threads 24 chr${i}.vcf.gz; done
# bcftools concat --threads 24 -a -D chr{1..22}.vcf.gz chrX.vcf.gz -Oz -o camgwas.vcf.gz
# 
# for i in bantu semibantu; do for j in {1..22}; do bcftools view -S ${i}.bcf.ids --force-samples chr${j}.vcf.gz -Oz -o ${i}_chr${j}.vcf.gz; done; done
