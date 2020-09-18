#!/usr/bin/env bash

sed -i 's/ - /_/g' malgenCP1_Bu_updated2.csv 
cut -f2,3,7,14-21 -d',' malgenCP1_Bu_updated2.csv > malgenCP1_Bu_updated_no_space_selected_phenos.csv
Rscript -e 'a <- read.table("malgenCP1_Bu_updated_no_space_selected_phenos.csv", h=T, fill=T, sep=","); f <- read.table("all_illum.tsv", h=T); fa <- merge(f,a,by="malID"); write.table(fa,"malgenCP1_Bu_updated_phenos.tsv",col.names=T, row.names=F, quote=F, sep="\t")'
