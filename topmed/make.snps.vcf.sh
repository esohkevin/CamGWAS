#!/usr/bin/env bash

for i in *.snps; do grep -v CHR $i | awk '{print $1"\t"$2"\t"$3"\t"$4}' | sed 's/\//\t/g' | sort -g | sed '1 i ##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT' | bgzip -c > ${i}.vcf.gz; done

for j in sbba sb ba; do for i in ${j}.*.snps; do grep -v CHR $i | awk '{print $1"\t"$2"\t"$3"\t"$4}' | sed 's/\//\t/g'; done | sort -g | uniq | sed '1 i ##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT' | bgzip -c > ${j}.snps.vcf.gz; done

