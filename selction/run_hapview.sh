#!/usr/bin/env bash

hv="${HOME}/bioTools/Haploview.jar"

for pop in sbba sb ba; do
    java -jar ${hv} -n -log logfile.txt -pedfile ${pop}.hv.ped -info ${pop}.hv.info -out ${pop}.out -chromosome 11 -assocCC -permtests 10000 -minMAF 0.01 -check -ldvalues ${pop}.ld -ldcolorscheme DEFAULT -infoTrack -png -memory 2000
done
