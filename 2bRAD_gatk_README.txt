Full-blown GATK pipeline for 2bRAD data, based on bowtie2 mapping to a reference genome
September 14, 2014, Mikhail Matz matz@utexas.edu

Includes analysis (making data for genome scans, running BayeScan and ADMIXTURE)

The idea is to copy the chunks separated by empty lines below and paste them into your cluster 
terminal window consecutively. 

The lines beginning with hash marks (#) are explanations and additional instructions - 
please make sure to read them before copy-pasting.
(some of these are for previously held workshops, ecogeno2013 and mega2014) 

In addition to the scripts coming with this distribution,
you will need the following software installed and available (note: TACC already has them as modules):

python: http://www.python.org/getit/
fastx_toolkit: http://hannonlab.cshl.edu/fastx_toolkit/download.html
bowtie2: http://bowtie-bio.sourceforge.net/index.shtml 
samtools: http://sourceforge.net/projects/samtools/files/
picard: http://sourceforge.net/projects/picard/files/
gatk: http://www.broadinstitute.org/gatk/download
vcftools: http://vcftools.sourceforge.net/ 
PGDspider2: http://www.cmpg.unibe.ch/software/PGDSpider/

NOTE: this walkthrough has been written for lonestar cluster of the Texas
Advanced Computer Center, which has 12 cores per node and uses Sun Grid Engine 
(SGE) batch processing system. To adopt the walkthough to your cluster you 
would need to edit the launcher_creator.py script to make its default settings
compatible with your cluster. 

==============================================

# concatenating the reads files (make sure to unzip them first!):
# NOTE: Only run this chunk if your samples are spread over several lanes.
# Assuming your filed have the extension fastq (edit that in the line below if not),
# replace "2b_(.+)_L00" below with the actual pattern to recognize sample
# identifier in the file name. The pattern above will work for file names such as
# Sample_2b_K208_L007_R1.cat.fastq (the identifier is K208)
ngs_concat.pl fastq "2b_(.+)_L00"


# Trimming the reads, splitting them by secondary barcode and removing PCR duplicates
# The sampleID parameter here, '3', defines a chunk in the filename (separated
# by underscores) that is to be kept; 
# for example, for a filename like this: Sample_2b_M11_L006_R1.cat.fastq 
# it would be reasonable to specify 'sampleID=3' to keep only 'M11' as a 
# sample identifier. If the sampleID parameter is not specified, the whole filename up to the first dot will be kept.
# NOTE: if you ran ngs_concat.pl (above), run not the next but the second-next line (after removing # symbol) :
2bRAD_trim_launch_dedup.pl fastq sampleID=3 > trims
launcher_creator.py -j trims -n trims -l trimjob
qsub trimjob

# NB: use this command instead of the one above if you have 2bRAD libraries without 
# degenerate tag but with 4-base in-line barcode:
2bRAD_trim_launch.pl fastq barcode2=4 > trims
# And if you have the original-style libraries without degenerate tag and without in-line barcode, use this:
2bRAD_trim_launch.pl fastq > trims


# quality filtering using fastx_toolkit
module load fastx_toolkit
ls *.tr0 | perl -pe 's/^(\S+)\.tr0$/cat $1\.tr0 \| fastq_quality_filter -q 20 -p 90 >$1\.trim/' >filt0

# NOTE: run the next line ONLY if your qualities are 33-based 
# (if you don't know just try to see if it works eventually, if you get errors from fastx_toolkit, try the other one):
	cat filt0 | perl -pe 's/filter /filter -Q33 /' > filt
#if you did NOT run the line above, run this one (after removing # symbol):
#	mv filt0 filt

launcher_creator.py -j filt -n filt -l filtjob
qsub filtjob

#----------------------------------
# setting up the genome reference

# concatenating genome contigs into small number of "pseudo-chromosomes" (to help memory usage in GATK)
concatFasta.pl fasta=mygenome.fasta

# normalize fasta records to ensure all lines are of the same length, using Picard
module load picard
java -Xmx1g -jar $TACC_PICARD_DIR/NormalizeFasta.jar INPUT=mygenome_cc.fasta OUTPUT=mygenome_ccn.fasta

# Now, edit these lines THROUGHOUT THE TEXT to fit the location and name of your genome reference 
export GENOME_FASTA=mygenome_ccn.fasta  # for mega2014 class, this would be amil_genome_fold_c_cc.fasta
export GENOME_DICT=mygenome_ccn.dict  # same name as genome fasta but with .dict extension; mega2014: amil_genome_fold_c.dict
export GENOME_PATH=where/genome/is/
export GENOME_REF=where/genome/is/mygenome_ccn.fasta

# ecogeno2014: skip this
# UNLESS you are working on TACC, edit these accordingly and execute:
export TACC_GATK_DIR=/where/gatk/is/installed/
export TACC_PICARD_DIR=/where/picard/is/installed/
# NOTE that you will have to execute the above two lines every time you re-login!


module load bowtie/2.1.0

# creating genome indexes:
cd $GENOME_PATH
echo 'bowtie2-build $GENOME_FASTA $GENOME_FASTA' >btb
launcher_creator.py -j btb -n btb -l btbl
qsub btbl

module load samtools
samtools faidx $GENOME_FASTA

module load picard
# edit this line so your .dict filename is the same as genome fasta filename
# only with .dict extension instead of .fasta
export GENOME_DICT=amil_genome_fold_c_cc.dict 
java -jar $TACC_PICARD_DIR/CreateSequenceDictionary.jar R=$GENOME_FASTA  O=$GENOME_DICT

cds
cd rad/gatk/

#----------------------------------
# mapping reads and reformatting mapped data files

# aligning with bowtie2 :
module load bowtie/2.1.0
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_REF > bt2
launcher_creator.py -j bt2 -n maps -l bt2.job -t 3:00:00 -q normal
# calculate the max possible number of cores reasonable: ceiling(Nsamples/wayness)*12  
cat bt2.job | perl -pe 's/12way \d+/4way 216/' > bt2l.job
qsub bt2l.job

ls *.bt2.sam > sams
cat sams | wc -l  
# do you have sams for all your samples?... If not, rerun the chunk above

# making bam files
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
module load picard
module load samtools
module load jdk64
cat sams | perl -pe 's/(\S+)\.sam/samtools import \$GENOME_REF $1\.sam $1\.unsorted\.bam && samtools sort $1\.unsorted\.bam $1\.sorted && java -Xmx5g -jar \$TACC_PICARD_DIR\/AddOrReplaceReadGroups\.jar INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b
launcher_creator.py -j s2b -n s2b -l s2b.job -q normal
cat s2b.job | perl -pe 's/12way \d+/4way 216/' > sam2bam.job # REPLACE NNN WITH 12*ceiling([number of sam files]/4)
qsub sam2bam.job

rm *sorted*
ls *bt2.bam > bams
cat bams | wc -l  
# do you have bams for all your samples?... If not, rerun the chunk above

#----------------------------------
# WARNING!!! DO NOT RUN THIS FOR RAD!!! 
# marking duplicate reads 
# module load picard
# export GENOME_REF=where/genome/is/mygenome_ccn.fasta
#cat bams | perl -pe 's/(\S+)\.bam/java -Xmx4g -jar \$TACC_PICARD_DIR\/MarkDuplicates\.jar INPUT=$1\.bam OUTPUT=$1\.dedup.bam METRICS_FILE=$1\.metrics AS=true CREATE_INDEX=true/' > dd
#launcher_creator.py -j dd -n dd -l dd.job -q normal
#cat dd.job | perl -pe 's/12way \d+/4way NNN/' > ddup.job  # REPLACE NNN WITH 12*ceiling([number of bam files]/4)
#qsub ddup.job
#ls *.dedup.bam > bams

#----------------------------------
# starting GATK
# realigning around indels:

# step one: finding places to realign:
module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx5g -jar \$TACC_GATK_DIR\/GenomeAnalysisTK\.jar -T RealignerTargetCreator -R \$GENOME_REF -I $1\.bam -o $1\.intervals/' >intervals
launcher_creator.py -j intervals -n intervals -l inter.job -q normal
cat inter.job | perl -pe 's/12way \d+/4way 216/' > intervals.job  # if you are not EcoGeno2014, replace 216 with 12*ceiling([number of bam files]/4)
qsub intervals.job

# did it run for all files? is the number of *.intervals files equal the number of *.bam files?
# if not, rerun the chunk above
ll *.intervals | wc -l

# step two: realigning
module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx5g -jar \$TACC_GATK_DIR\/GenomeAnalysisTK\.jar -T IndelRealigner -R \$GENOME_REF -targetIntervals $1\.intervals -I $1\.bam -o $1\.real.bam -LOD 0\.4/' >realign
launcher_creator.py -j realign -n realign -l realig.job -q normal
cat realig.job | perl -pe 's/12way \d+/4way 216/' > realign.job # if you are not EcoGeno2014, replace 216 with 12*ceiling([number of bam files]/4)
qsub realign.job

# did it run for all files? is the number of *.intervals files equal the number of *.bam files?
# if not, rerun the chunk above
ll *.real.bam | wc -l

#----------------------------------
# launching GATK UnifiedGenotyper for round 1 (about 30 min)
# note: it is a preliminary run needed for base quality recalibration,
module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
ls *real.bam > bams
echo '#!/bin/bash
#$ -V
#$ -cwd
#$ -N unig
#$ -A mega2014
#$ -pe 1way 12
#$ -q normal   
#$ -l h_rt=24:00:00
#$ -M matz@utexas.edu
#$ -m be
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R $GENOME_REF -nt 12 -nct 1 \' >unig
cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unig
echo '-o round1.vcf ' >> unig
qsub unig

#----------------------------------
# base quality score recalibration (BQSR)

# creating high-confidence (>75 quality percentile) snp sets for 
# base quality recalibration using Kyle Hernandez's tool. (ignore the warnings)
GetHighQualVcfs.py  -i round1.vcf --percentile 75 -o .

# recalibrating quality scores
# step one: creating recalibration reports
module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
ls *real.bam > bams
cat bams | perl -pe 's/(\S+)\.real\.bam/java -Xmx20g -jar \$TACC_GATK_DIR\/GenomeAnalysisTK\.jar -T BaseRecalibrator -R \$GENOME_REF -knownSites $1_HQ\.vcf -I $1\.real\.bam -o $1\.real\.recalibration_report.grp/' >bqsr
launcher_creator.py -j bqsr -n bqsrHC -l bqsrjob -q normal
cat bqsrjob | perl -pe 's/12way \d+/1way 828/' > bqsr.job # if you are not EcoGeno2014, replace 828 with 12*[number of bam files]
qsub bqsr.job

# did it run for all files? is the number of *.grp files equal the number of *.real.bam files?
# if not, rerun the chunk above
ll *.real.bam | wc -l
ll *.grp | wc -l

# step two: rewriting bams according to recalibration reports
module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
cat bams | perl -pe 's/(\S+)\.bam/java -Xmx10g -jar \$TACC_GATK_DIR\/GenomeAnalysisTK\.jar -T PrintReads -R \$GENOME_REF -I $1\.bam -BQSR $1\.recalibration_report.grp -o $1\.recal\.bam /' >bqsr2
launcher_creator.py -j bqsr2 -n bqsr2 -l bqsrjob2 -q normal
cat bqsrjob2 | perl -pe 's/12way \d+/2way 420/' > bqsr2.job # if you are not EcoGeno2014, replace 420 with 12*ceiling([number of bam files]/2)
qsub bqsr2.job

# did it run for all files? is the number of *.recal.bam files equal the number of *.real.bam files?
# if not, rerun the chunk above
ll *.real.bam | wc -l
ll *.recal.bam | wc -l

ls *.recal.bam > bams

#----------------------------------
# Second iteration of UnifiedGenotyper (on quality-recalibrated files)
# this time FOR REAL! 
# in you need indels, run the same process separately with --genotype_likelihoods_model INDEL
# I do not recommend indel tracing for 2bRAD since the tags are too short for confident indels. 
# If you still want to try, note that the subsequent recalibration stages would 
# have do be done separately for indels, 

module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
ls *.recal.bam > bams
echo '#!/bin/bash
#$ -V
#$ -cwd
#$ -N unig2
#$ -A mega2014
#$ -pe 1way 12
#$ -q normal   
#$ -l h_rt=24:00:00
#$ -M matz@utexas.edu
#$ -m be
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R $GENOME_REF -nt 12 -nct 1 \
--genotype_likelihoods_model SNP \' >unig2
cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unig2
echo '-o round2.vcf ' >> unig2
qsub unig2

# renaming samples in the vcf file, to get rid of trim-shmim etc
cat round2.vcf | perl -pe 's/\.fq\.trim\.bt2//g' | perl -pe 's/\.trim\.bt2//g' | perl -pe 's/^chrom/chr/' >round2.names.vcf

#----------------------------------
# variant quality score recalibration (VQSR)

# extracting SNPs that are consistently genotyped in replicates, polymorphic, 
# and have the fraction of alternative allele not too low and not too close to 0.5:
replicatesMatch.pl vcf=round2.names.vcf replicates=clonepairs.tab polyonly=1 >vqsr.vcf

# determining transition-transversion ratio for true snps (will need it for tranche calibration)
vcftools --vcf vqsr.vcf --TsTv-summary
# Ts/Tv ratio: 1.455  # put your actual number into the next code chunk, --target_titv

# Recalibrating genotype calls: VQSR
# step one - creating recalibration models (30 sec)
module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T VariantRecalibrator \
-R $GENOME_REF -input round2.names.vcf -nt 12 \
-resource:repmatch,known=true,training=true,truth=true,prior=10  vqsr.vcf \
-an QD -an MQ -an FS -mode SNP --maxGaussians 4 \
--target_titv 1.455 -tranche 85.0 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile round2.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R 

# fixing the R script (assumes outdated version of ggplot2):
cat recalibrate_SNP_plots.R | perl -pe 's/opts\(/theme\(/'g | perl -pe 's/theme_/element_/g' | perl -pe 's/\+ theme\(title=\"model PDF\"\)//g'  >recalibrateSNPs.R
# now copy al recalibrate* files to your laptop, run the R script, examine the resulting plot and tranches.pdf

# applying recalibration:
module load gatk/3.1.1
module load jdk64
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $GENOME_REF -input round2.names.vcf -nt 12 \
--ts_filter_level 95.0 -mode SNP \
-recalFile round2.recal -tranchesFile recalibrate_SNP.tranches -o gatk_after_vqsr.vcf

#----------------------------------
# ALTERNATIVE to GATK-based recalibration (if tranches are weird and error models fail to converge):
# non-parametric quantile-based recalibration a-la de novo pipeline

replicatesMatch.pl vcf=round2.names.vcf replicates=clonepairs.tab polyonly=1 >vqsr.vcf
recalibrateSNPs_gatk.pl vcf=round2.names.vcf true=vqsr.vcf >gatk_after_vqsr.vcf

# your best quality filter setting is the one with maximum "gain"

#----------------------------------

# restoring original contig names and coordinates:
echo 'retabvcf.pl vcf=gatk_after_vqsr.vcf tab=where/genome/is/mygenome_cc.tab > retab.vcf '>rt
launcher_creator.py -j rt -n rt -l rtj
qsub rtj

# thinning: leaving one SNP per tag, the one with max minor allele frequency
thinner.pl vcf=retab.vcf >thin.vcf

# applying filter and selecting polymorphic biallelic loci genotyped in 80% or more individuals
# for parametric (GATK-based) recalibration"
vcftools --vcf thin.vcf --remove-filtered-all --max-missing 0.8  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filt0
#After filtering, kept 18792 out of a possible 95211 Sites

# for non parametric (GATK-based) recalibration: replace --minQ 15 in the following line
# with the quantile of the highest "gain" as reported by recalibrateSNPs_gatk.pl
vcftools --vcf thin.vcf --minQ 15 --max-missing 0.8  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filt0

# genotypic match between pairs of replicates 
# (the most important one is the last one, Heterozygote Discovery Rate)	
repMatchStats.pl vcf=filt0.recode.vcf replicates=clonepairs.tab 
#pair	gtyped	match	[ 00	01	11 ]	HetMatch	HomoHetMismatch	HetNoCall	HetsDiscoveryRate
#K210:K212	18348	18116(98.7%)	 [70%	17%	13% ]	3130	194	5	0.97	
#K212:K213	18433	18149(98.5%)	 [70%	17%	13% ]	3145	200	3	0.97	
#K213:K216	18414	17950(97.5%)	 [70%	17%	13% ]	3052	276	8	0.96	
#K211:K219	18401	18176(98.8%)	 [69%	18%	13% ]	3321	154	1	0.98	

# creating final filtered file without clones (must list them in the file clones2remove):
# (for non-parametric recalibration, replace --remove-filtered-all with --minQ [quantile of the highest gain] )
vcftools --vcf thin.vcf --remove clones2remove --remove-filtered-all --max-missing 0.8  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out OK
# After filtering, kept 19285 out of a possible 95211 Sites

# creating unthinned dataset for Tajima's D calculations
# (for non-parametric recalibration, replace --remove-filtered-all with --minQ [quantile of the highest gain] )
vcftools --vcf retab.vcf --remove clones2remove --remove-filtered-all --max-missing 0.8  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out OKunthin
# After filtering, kept 36677 out of a possible 174106 Sites

# ===================================================
# Analysis

# creating population lists
# in the example below all the files beginning with O are from one pop, K from another, so
ls K*.recal.bam | perl -pe 's/\..+//' >K.pop
ls O*.recal.bam | perl -pe 's/\..+//' >O.pop
cat O.pop K.pop  | perl -pe 's/((.).+)/$1\t$2/' >OK.pops

#-------------
# need PGDspider - the format-converting software!
# install PDGspider on your laptop: 
# http://www.cmpg.unibe.ch/software/PGDSpider/#Download_and_Installation_Instructions  

# install PGDspider on lonestar:
cd ~/bin
wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.1.zip
unzip PGDSpider_2.0.7.1.zip
cd -

#-------------
# BayeScan: searching for Fst outliers

# converting vcf to bayescan format

# creating configuration file for PDGspider
nano vcf2bayescan.spid

# paste this (edit if needed):
####################################
# VCF Parser questions
PARSER_FORMAT=VCF
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=./OK.pops
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=

# GESTE / BayeScan Writer questions
WRITER_FORMAT=GESTE_BAYE_SCAN
# Specify which data type should be included in the GESTE / BayeScan file  (GESTE / BayeScan can only analyze one data type per file):
GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP
####################################

java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar -inputfile OK.recode.vcf -outputfile OK.bayescan -spid vcf2bayescan.spid 

# Download and install BayeScan 
cd ~/bin
wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip
cp BayeScan2.1/binaries/BayeScan2.1_linux64bits bayescan
chmod +x bayescan
rm -r BayeScan*

echo 'bayescan OK.bayescan -threads=192' >bs
launcher_creator.py -j bs -n bs -l bsl
cat bsl | perl -pe 's/12way 12/12way 192/' | perl -pe 's/h_rt=1/h_rt=24/' |  perl -pe 's/development/normal/' >bsll
qsub bsll
# this one takes a few hours

#-------------
# pi: (nucleotide diversity = expected heterozygosity)

vcftools --vcf OK.recode.vcf --keep O.pop --site-pi
mv out.sites.pi O.pi
vcftools --vcf OK.recode.vcf  --keep K.pop --site-pi
mv out.sites.pi K.pi

#-------------
# Fst:

vcftools --vcf OK.recode.vcf  --weir-fst-pop O.pop --weir-fst-pop K.pop  
# 0.011437
mv out.weir.fst OK.fst

#-------------
# Tajima's D (separately for two pops):

vcftools --vcf OKunthin.recode.vcf  --keep K.pop --TajimaD 75000
mv out.Tajima.D K.td
vcftools --vcf OKunthin.recode.vcf  --keep O.pop --TajimaD 75000
mv out.Tajima.D O.td

#-------------
# LD (using only bi-allelic sites within 15 kb of each other) - these two jobs take quite some time...
echo 'vcftools --vcf OK.recode.vcf --geno-r2 --ld-window-bp 15000' >ld
launcher_creator.py -j ld -n ld -l ldj
cat ldj | perl -pe 's/12way/1way/' | perl -pe 's/h_rt=1/h_rt=12/' |perl -pe 's/development/normal/' > ldjj
qsub ldjj

# wait for it to complete (~15 min), then do this (another 15-20 min):
module load R
echo 'Rscript ~/bin/analyzeLD.R' > ldr
launcher_creator.py -j ldr -n ldr -l ldrj
cat ldrj | perl -pe 's/12way/1way/' | perl -pe 's/h_rt=1/h_rt=24/' |perl -pe 's/development/normal/' > ldrjj
qsub ldrjj

# scp (or WinSCP) the resulting files: OK.baye_fst.txt, O.pi, K.pi, OK.fst, K.td, O.td, LD_quantized.RData 
# to your laptop

# use genomeScanPlots.R to plot it

#-------------
# ADMIXTURE
#-------------

# installing ADMIXTURE
cd ~/bin/
wget https://www.genetics.ucla.edu/software/admixture/binaries/admixture_linux-1.23.tar.gz --no-check-certificate
tar vxf admixture_linux-1.23.tar.gz 
mv admixture_linux-1.23/admixture .
cd -

# installing plink 1.09:
cd ~/bin
wget https://www.cog-genomics.org/static/bin/plink140918/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip
cd-

# for denovo vcf, creating a dataset with fake "chr" chromosome designations
cat denovo.vcf | perl -pe 's/locus(\d)(\d+)\t\d+/chr$1\t$2/' >chrom.recode.vcf

# for GATK-based vcf, creating a thinned but not re-tabbed dataset with "chr" chromosome designations
thinner.pl vcf=gatk_after_vqsr.vcf | perl -pe 's/chrom/chr/g' >thinchrom.vcf
vcftools --vcf thinchrom.vcf --remove clones2remove --remove-filtered-all --max-missing 0.8  --min-alleles 2 --max-alleles 2 --recode --out chrom

# reformatting VCF into plink binary BED format
plink --vcf chrom.recode.vcf --make-bed --out admix

# running ADMIXTURE with
# cross-validation to select K (bash script)
for K in 1 2 3 4 5 6; \
do admixture --cv admix.bed $K | tee log${K}.out; done

grep -h CV log*.out # select K with the minimal value

# listing individuals from vcf file
grep "#CHROM" chrom.recode.vcf | perl -pe 's/\t/\n/g' | grep [0-9] > inds.list

# scp the *.Q and inds.list files to laptop, plot it in R:
# use admixturePlotting.R to plot (will require minor editing - population names)
