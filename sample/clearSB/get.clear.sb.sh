#!/usr/bin/env bash

for i in $(cut -f1 -d' ' ../../../samples/include_semibantu.txt); do grep -w $i ../malgenCP1_Bu_updated_phenos.tsv; done | awk '$8!="NA" && $4=="CASE"' | sed '1 i malID\tgtc_id\tsampleID.x\tcaseorcontrol\testimated_sex\tcurated_ethnicity\tsampleID.y\tage_in_months\tparacat\tCC\tSMA\thyperpyrexia\thyperparasitaemia\tSM\tUM\tage_under5' > sb.clear43cases.txt

for i in $(cut -f1 -d' ' ../../../samples/include_semibantu.txt); do grep -w $i ../malgenCP1_Bu_updated_phenos.tsv; done | awk '$8!="NA" && $4!="CASE"' | sed '1 i malID\tgtc_id\tsampleID.x\tcaseorcontrol\testimated_sex\tcurated_ethnicity\tsampleID.y\tage_in_months\tparacat\tCC\tSMA\thyperpyrexia\thyperparasitaemia\tSM\tUM\tage_under5' > sb.clear115controls.txt


grep -w "SEMI_BANTU" sb.clear43cases.txt | awk '{print $2,$2}' | shuf -n 50 -o sb.clear50cases.txt
grep -w "SEMI_BANTU" sb.clear115controls.txt | awk '{print $2,$2}' | shuf -n 50 -o sb.clear50controls.txt
