#!/usr/bin/env bash

# plink2 --bfile chr10 --export oxford bgz --chr 10 --out chr10.ox --covar ../../analysis/assoc/sbba.pca.glm.cov 
# 
# cut -f5- chr10.ox.cov | sed 's/\t/ /g' | sed '2 i C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C' > ox.cov
# 
# paste chr10.ox.sample ox.cov | sed 's/\t/ /g' > chr10.ox.cov
# 
# mv chr10.ox.cov chr10.ox.sample
# 
# rm ox.cov 

snptest_v2.5.4-beta3 -data chr10.ox.gen.gz chr10.ox.sample -o chr10.snptest.assoc -bayesian 1 -pheno PHENO1 -method score -filetype gen
