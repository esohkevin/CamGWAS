###### CAMEROON GWAS DATA ANALYSIS PIPELINE ######
# Data: 
#	Genotype:	camgwas_merge.vcf.gz
#			${raw}-1.gen 
#
#			from: 

#runplink() {
images="../images/"
samples="../samples/"
mkdir -p ../images
king="${HOME}/bin/king"
raw['test']=NULL

#-- Check for duplicate SNPs
plink \
	--bfile ${raw} \
	--allow-no-sex \
	--list-duplicate-vars ids-only suppress-first \
	--out dups

#-- by the PARs using the b37 coordinates while removing duplicate SNPs
plink \
	--bfile ${raw} \
	--make-bed \
	--exclude dups.dupvar \
	--split-x b37 \
	--allow-no-sex \
	--out temp1

#-- Update SNPID names with rsids
cut -f2 temp1.bim > all.snps.ids
cut -f1 -d',' all.snps.ids > all.rs.ids
paste all.rs.ids all.snps.ids > allMysnps.txt

plink \
        --bfile temp1 \
        --update-name allMysnps.txt 1 2 \
        --allow-no-sex \
        --make-bed \
        --out ${raw}-1

#-- LD-prune the raw data before sex check
plink \
        --bfile ${raw}-1 \
        --allow-no-sex \
        --indep-pairwise 50 10 0.2 \
	--set-hh-missing \
        --out prunedsnplist

#-- Now extract the pruned SNPs to perform check-sex on
plink \
        --bfile ${raw}-1 \
        --allow-no-sex \
        --extract prunedsnplist.prune.in \
        --make-bed \
        --out check-sex-data

#-- Check for sex concordance
plink \
	--bfile check-sex-data \
	--check-sex \
	--set-hh-missing \
	--allow-no-sex \
	--out check-sex-data

#-- Extract FIDs and IIDs of individuals flagged with error 
#-- (PROBLEM) in the .sexcheck file (failed sex check)
grep "PROBLEM" check-sex-data.sexcheck | awk '{print $1"\t"$2}' > fail-checksex.qc

#-- Compute missing data stats
plink \
	--bfile ${raw}-1 \
	--missing \
	--allow-no-sex \
	--set-hh-missing \
	--out ${raw}-1

#-- Compute heterozygosity stats
plink \
	--bfile ${raw}-1 \
	--het \
	--allow-no-sex \
	--set-hh-missing \
	--out ${raw}-1

echo -e """\e[38;5;40m
	##########################################################################
	##	    Perform per individual missing rate QC in R			##
	##########################################################################
	\e[0m
	"""
echo -e "\n\e[38;5;40mNow generating plots for per individual missingness in R. Please wait...\e[0m"

Rscript indmissing.R ${raw}-1.het ${raw}-1.imiss

#-- Extract a subset of ${raw}-rel individuals to produce an IBD 
#-- report to check duplicate or related individuals baseDird on autosomes
plink2 \
	--bfile ${raw}-1 \
	--autosome \
	--maf 0.30 \
	--geno 0.02 \
	--hwe 1e-50 keep-fewhet \
	--allow-no-sex \
	--make-bed \
	--out ${raw}-rel

#-- Run KING software for relatedness
      $king \
      	-b ${raw}-rel.bed \
      	--ibdseg \
      	--rplot \
      	--prefix ${raw}-rel
      Rscript ${raw}-rel_ibd1vsibd2.R ${raw}-rel.segments.gz ibdseg
       awk '$10!="3rd" && $10!="4th" && $10!="UN"' ${raw}-rel.seg | cut -f1 | sort | uniq | awk '{print $1"\t"$1}' > ${raw}-unrelated.ids

#-- Prune the list of ${raw}-rel SNPs to remove those that fall within 
#-- 50bp with r^2 > 0.2 using a window size of 5bp
plink \
	--bfile ${raw}-rel \
	--allow-no-sex \
	--indep-pairwise 50 10 0.2 \
	--out prunedsnplist

#-- Now generate the IBD report with the set of pruned SNPs 
#-- (prunedsnplist.prune.in - IN because they're the ones we're interested in)
plink \
	--bfile ${raw}-rel \
	--allow-no-sex \
	--extract prunedsnplist.prune.in \
	--genome \
	--out caseconpruned

echo -e """\e[38;5;40m
	#########################################################################
	#              Perform IBD analysis (relatedness) in R                  #
	#########################################################################
	\e[0m
	"""
echo -e "\n\e[38;5;40mNow generating plots for IBD analysis in R. Please wait...\e[0m"

#Rscript ibdana.R

#-- Merge IDs of all individuals that failed per individual qc
#cat fail-checksex.qc  fail-het.qc  fail-mis.qc duplicate.ids1 | sort | uniq > fail-ind.qc
cat fail-checksex.qc  fail-het.qc  fail-mis.qc ${raw}-unrelated.ids | sort | uniq > fail-ind.qc

#-- Remove individuals who failed per individual QC
plink \
	--bfile ${raw}-1 \
	--make-bed \
	--allow-no-sex \
	--set-hh-missing \
	--remove fail-ind.qc \
	--out ${raw}-ind-qc

#-- Per SNP QC
#-- Compute missing data rate for ${raw}-ind-qc data
plink \
	--bfile ${raw}-ind-qc \
	--allow-no-sex \
	--set-hh-missing \
	--missing \
	--out ${raw}-ind-qc

#-- Compute MAF
plink \
	--bfile ${raw}-ind-qc \
	--allow-no-sex \
	--set-hh-missing \
	--freq \
	--out ${raw}-ind-qc

#-- Compute differential missing genotype call rates (in cases and controls)
plink \
	--bfile ${raw}-ind-qc \
	--allow-no-sex \
	--set-hh-missing \
	--test-missing \
	--out ${raw}-ind-qc

echo -e """\e[38;5;40m
	#########################################################################
	#                        Perform per SNP QC in R                        #
	#########################################################################
	\e[0m
	"""
echo -e "\n\e[38;5;40mNow generating plots for per SNP QC in R. Please wait...\e[0m"

Rscript snpmissing.R ${raw}-ind-qc.lmiss ${raw}-ind-qc.frq ${raw}-ind-qc.missing

#-- Remove SNPs that failed per marker QC
plink2 \
	--bfile ${raw}-ind-qc \
	--exclude fail-diffmiss.qc \
	--allow-no-sex \
	--maf 0.0001 \
	--hwe 1e-50 keep-fewhet \
	--geno 0.05 \
	--make-bed \
	--merge-x \
	--out qc-${raw}

echo -e """\e[38;5;40m
	#########################################################################
	#                          ChrX Quality Control                         #
	#########################################################################
	\e[0m
	"""
echo -e "\n\e[38;5;40mNow generating plots for per SNP QC in R. Please wait...\e[0m"

#-- Extract only autosomes for subsequently merging with QCed chrX
plink \
        --bfile qc-${raw} \
        --allow-no-sex \
        --make-bed \
        --out qc-${raw}

plink \
	--bfile qc-${raw} \
	--allow-no-sex \
        --set-hh-missing \
	--make-bed \
	--autosome \
	--out qc-${raw}-autosome

#-- Extract only chrX for QC
plink \
	--bfile qc-${raw} \
	--allow-no-sex \
        --set-hh-missing \
	--make-bed \
	--chr 23 \
	--out qc-${raw}-chrX

#-- Compute differential missingness
plink \
        --bfile qc-${raw}-chrX \
        --allow-no-sex \
        --set-hh-missing \
        --test-missing \
        --out qc-${raw}-chrX

echo -e """\e[38;5;40m
	#########################################################################
	#                          chrX per SNP QC in R                         #
	#########################################################################
	\e[0m
	"""
echo -e "\n\e[38;5;40mPerforming ChrX per SNP QC in R. Please wait...\e[0m"

#Rscript xsnpmissing.R
awk '$5<1e-8' qc-${raw}-chrX.missing > fail-Xdiffmiss.qc

#-- Now remove SNPs that failed chrX QC
plink2 \
        --bfile qc-${raw}-chrX \
        --exclude fail-Xdiffmiss.qc \
        --allow-no-sex \
        --maf 0.0001 \
        --hwe 1e-50 keep-fewhet \
        --geno 0.05 \
        --make-bed \
        --max-alleles 2 \
	--out qc-${raw}-chr23 

#-- Merge autosome and chrX data sets again
plink \
	--bfile qc-${raw}-chr23 \
	--make-bed \
	--set-hh-missing \
	--out qc-${raw}-chr23

plink \
	--bfile qc-${raw}-autosome \
	--allow-no-sex \
	--bmerge qc-${raw}-chr23 \
	--set-hh-missing \
	--out qc-${raw}
#done
rm *~

echo -e """\e[38;5;40m
	#########################################################################
	#                     	   Updating QC rsids                            #
	#########################################################################
	\e[0m
	"""
cut -f1,4 qc-${raw}.bim | \
	sed 's/\t/:/g' > qc-${raw}.pos
cut -f2 qc-${raw}.bim > qc-${raw}.ids
paste qc-${raw}.ids qc-${raw}.pos > qc-${raw}-ids-pos.txt

plink \
	--bfile qc-${raw} \
	--update-name qc-${raw}-ids-pos.txt 2 1 \
	--allow-no-sex \
	--make-bed \
	--set-hh-missing \
	--exclude-snp kgp21103953 \
	--out qc-${raw}

plink \
	--bfile qc-${raw} \
	--update-name ../../../db/updateName.txt 1 2 \
	--allow-no-sex \
	--make-bed \
	--set-hh-missing \
	--out qc-${raw}

echo -e """\e[38;5;40m
	#########################################################################
	#                     Run Imputation Prep Script                        #
	#########################################################################
	\e[0m
	"""

