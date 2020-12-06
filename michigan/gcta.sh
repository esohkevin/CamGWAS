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
gc="${an}gcta/"
gcta="${HOME}/bin/gcta64"


${gcta} --bfile cam.clean --threads 24  --ld-score-region 100 --out cam

./gcta.r

seq 4 | parallel echo "--bfile cam.clean --extract snp_group{}.txt --threads 24 --make-grm --out cam_group{}" | xargs -I input -P5 sh -c "${gcta} input"
#gcta64 --bfile cam.clean --extract snp_group2.txt --make-grm --out cam_group2
#gcta64 --bfile cam.clean --extract snp_group1.txt --make-grm --out cam_group3
#gcta64 --bfile cam.clean --extract snp_group2.txt --make-grm --out cam_group4

echo -e "cam_group1\ncam_group2\ncam_group3\ncam_group4\n" > multi_GRMs.txt

${gcta} --reml --mgrm multi_GRMs.txt --pheno ${gc}sbba.pca.grm.pheno --threads 24 --qcovar ${gc}sbba.pca.grm.cov --covar ${gc}sbba.cov.sex --reml-maxit 10000  --out cam

#--reml-no-constrain
