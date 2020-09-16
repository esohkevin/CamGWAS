#!/usr/bin/env bash

for pop in sbba sb ba; do 
	plink \
		--vcf hbb.vcf.gz \
		--recode HV-1chr \
		--keep-allele-order \
		--out ${pop}.hv \
		--chr 11 \
		--from-kb 5200 \
		--to-kb 5300 \
		--update-sex raw-camgwas.fam 3 \
		--pheno raw-camgwas.fam \
		--mpheno 4 \
		--double-id \
		--keep ${pop}.sample.ids
done
