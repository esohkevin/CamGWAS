#!/usr/bin/env bash

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


for p in sbba sb ba; do echo $p; done | parallel echo "--bfile cam.clean --recode 12 transpose --keep ${asd}{}.pca.glm.cov --threads 24 --keep-allele-order --maf 0.01 --autosome --out {}" | xargs -I input -P5 sh -c "plink input"
#  #for j in *.tfam; do echo "<<< ${j/.tfam/.cov} >>>"; for i in $(awk '{print $1}' ${j}); do grep -w $i ${asd}pca.glm.cov; done > ${j/.tfam/.cov}; done
#  
#  cut -f1-4,6-25 -d' ' ${asd}sbba.pca.glm.cov | sed '1d' > sbba.cov
#  cut -f1-4,6-25 -d' ' ${asd}sb.pca.glm.cov | sed '1d' > sb.cov
#  cut -f1-4,6-25 -d' ' ${asd}ba.pca.glm.cov | sed '1d' > ba.cov
#  
for p in sbba sb ba; do echo $p; done | parallel echo 1 {} | xargs -n2 -P10 ${as}run_emmax.sh
for p in sbba sb ba; do for i in clean sm sma; do mv ${p}.${i}.ps ${p}.${i}.nocov.ps; done; done
for p in sbba sb ba; do echo ${p}; done | parallel echo 2 {} | xargs -n2 -P10 ${as}run_emmax.sh
for i in *.ps; do sed -i '1 i SNP\tBETA\tSE\tP' ${i}; done
for i in *.ps; do echo $i; done | parallel echo ${an}fdr.R {} | xargs -I input -P10 sh -c "Rscript input"
for i in *.ps.adj.txt; do mv $i ${i/.adj.txt/}; done
for i in *.ps; do echo $i; done | parallel echo ${an}emmax2plink_assoc.R {} ${asd}qc-cam.bim | xargs -I input -P10 sh -c "Rscript input"
for i in *.assoc.ps; do echo $i; done | parallel echo ${an}assoc.R {} | xargs -I input -P 10 sh -c "Rscript input"
rm *.tped *.tfam


