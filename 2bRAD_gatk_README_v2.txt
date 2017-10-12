2bRAD reference-based walkthrough, version 2
June 2017 Mikhail Matz (matz@utexas.edu)

Major differences from previous version:
	- indels realignment and base quality score recalibration is  removed (does not make a difference for 2bRAD)
	- not using launcher_creator.py script

The idea is to copy the chunks separated by empty lines below and paste them into your cluster 
terminal window consecutively. 

The lines beginning with hash marks (#) are explanations and additional instructions - 
please make sure to read them before copy-pasting. 

In addition to the scripts coming with this distribution,
you will need the following software installed and available 
(note: TACC already has them as modules):

python: http://www.python.org/getit/
fastx_toolkit: http://hannonlab.cshl.edu/fastx_toolkit/download.html
cd-hit: https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
vcftools: http://vcftools.sourceforge.net/ 

Occasionally, an example screen output is listed after the program calls below.

NOTE: this walkthrough has been written for a SLURM system on TACC (Texas Advanced 
Computer System). You might wish to look through suggested SLURM job headers and 
adjust them according to your cluster specs (such as, max number of processes per node, 
queue types, allocation budgeting)

==============================================
IMPORTANT: DO NOT SEQUENCE 2BRAD LIBRARIES ALONE ON A HISEQ4000 LANE! INVARIANT BASES
(ADAPTOR, RESTRICTION SITE) WILL NOT BE READ WELL. MIX IN 20% OF PHIX SAMPLES TO 
DIVERSIFY; OR MIX YOUR 2BRAD LIBRARIES WITH SOMEONE ELSE'S NOT-2BRAD LIBRARIES 
(JUST MAKE SURE THE BARCODES DO NOT OVERLAP) 

If you did and you have used BcgI enzyme, use 2bRAD_trim_launch_dedup2.pl for trimming
==============================================

# Step 1: Trimming, deduplicating and quality-filtering the reads

# how do the reads look before trimming (displaying top 25 DNA sequences in O9.fq):
head -50 O9.fq | grep -E "^[ATGCN]+$"

# creating a file of commands to run (assuming reads are in fastq files, one file 
# per sample. Combine read files from multiple lanes for the same sample using ngs_concat.pl)
2bRAD_trim_launch_dedup.pl fastq > trims

# NB: use this command instead of the one above if you have 2bRAD libraries without 
# degenerate tag but with 4-base in-line barcode:
2bRAD_trim_launch.pl fastq barcode2=4 > trims
# And if you have the original-style libraries without degenerate tag and without in-line barcode, use this:
2bRAD_trim_launch.pl fastq > trims
# AND IF you sequenced your double-barcoded, deduplicatable libraries on HiSeq 4000 alone
# (resulting in poor quality at restriction site and adaptor bases) and you have used BcgI enzyme, use this:
2bRAD_trim_launch_dedup2.pl fastq > trims

# next chunk writes down a "SLURM header" for this job. Make sure you adjust it to your 
# system (processes per node) and number of *.fastq files
#  - N is number of nodes, - n is number of processes to run (i.e., number of fastq files)
# - A is your computer budget allocation
echo '#!/bin/bash
#SBATCH -J trims
#SBATCH -n 68
#SBATCH -N 2
#SBATCH -p normal
#SBATCH -o trims.o%j
#SBATCH -e trims.e%j
#SBATCH -t 1:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
'> heads
# assembling a launching the job:
cat heads trims >trims.slurm
sbatch trims.slurm

# do we have expected number of *.tr0 files created?
ls -l *.tr0 | wc -l

# quality filtering using fastx_toolkit
module load fastx_toolkit
ls *.tr0 | perl -pe 's/^(\S+)\.tr0$/cat $1\.tr0 \| fastq_quality_filter -q 20 -p 90 >$1\.trim/' >filt0

# NOTE: run the next line ONLY if your qualities are 33-based 
# (if you don't know just try to see if it works eventually, if you get errors from fastx_toolkit, try the other one):
	cat filt0 | perl -pe 's/filter /filter -Q33 /' > filt
#if you did NOT run the line above, run this one (after removing # symbol):
#	mv filt0 filt

# next chunk writes down a "SLURM header" for this job. Make sure you adjust it to your 
# system (processes per node) and number of *.tr0 files
echo '#!/bin/bash
#SBATCH -J filt
#SBATCH -n 68
#SBATCH -N 2
#SBATCH -p normal
#SBATCH -o filt.o%j
#SBATCH -e filt.e%j
#SBATCH -t 1:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
'> heads
# assembling a launching the job:
cat heads filt >filt.slurm
sbatch filt.slurm

# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

#====================
# Edit these lines THROUGHOUT THE TEXT to fit the location and name of your genome reference 
export GENOME_FASTA=mygenome_ccn.fasta  
export GENOME_DICT=mygenome_ccn.dict  # same name as genome fasta but with .dict extension
export GENOME_PATH=where/genome/is/
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
#=====================

# setting up genome reference

# concatenating genome contigs into small number of "pseudo-chromosomes" (to help memory usage in GATK)
concatFasta.pl fasta=mygenome.fasta

# normalize fasta records to ensure all lines are of the same length, using Picard
module load picard
java -Xmx1g -jar $TACC_PICARD_DIR/NormalizeFasta.jar INPUT=mygenome_cc.fasta OUTPUT=mygenome_ccn.fasta

# UNLESS you are working on TACC, edit these accordingly and execute:
export TACC_GATK_DIR=/where/gatk/is/installed/
export TACC_PICARD_DIR=/where/picard/is/installed/
# NOTE that you will have to execute the above two lines every time you re-login!

module load bowtie

# creating genome indexes:
export GENOME_FASTA=mygenome_ccn.fasta  
export GENOME_PATH=where/genome/is/
cd $GENOME_PATH
echo 'bowtie2-build $GENOME_FASTA $GENOME_FASTA' >btb
# next chunk writes down a "SLURM header" for this job. 
echo '#!/bin/bash
#SBATCH -J btb
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p development
#SBATCH -o btb.o%j
#SBATCH -e btb.e%j
#SBATCH -t 1:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
'> heads
# assembling a launching the job:
cat heads btb >btb.slurm
sbatch btb.slurm

module load samtools
samtools faidx $GENOME_FASTA

module load picard
# edit this line so your .dict filename is the same as genome fasta filename
# only with .dict extension instead of .fasta
export GENOME_FASTA=mygenome_ccn.fasta  
export GENOME_DICT=mygenome_ccn.dict  
java -jar $TACC_PICARD_DIR/CreateSequenceDictionary.jar R=$GENOME_FASTA  O=$GENOME_DICT

#-------------
# Mapping and formatting bam files 

# map with bowtie2 with end-to-end matching

export GENOME_REF=where/genome/is/mygenome_ccn.fasta
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_REF | perl -pe 's/--score-min L,16,1 --local -L 16 //'> bt2
# adjust this SLURM header as needed
echo '#!/bin/bash
#SBATCH -J maps
#SBATCH -n 68
#SBATCH -N 2
#SBATCH -p normal
#SBATCH -o maps.o%j
#SBATCH -e maps.e%j
#SBATCH -t 4:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
'> heads
# assembling a launching the job:
cat heads bt2 >maps.slurm
sbatch maps.slurm

# what are mapping efficiencies? 
cat maps.e*

ls *.bt2.sam > sams
cat sams | wc -l  # number should match number of trim files

# making and sorting bam files (compressed sam files)

export GENOME_REF=where/genome/is/mygenome_ccn.fasta
module load picard-tools
module load samtools

cat sams | perl -pe 's/(\S+)\.sam/samtools import \$GENOME_REF $1\.sam $1\.unsorted\.bam && samtools sort -o $1\.sorted\.bam $1\.unsorted\.bam && java -Xmx5g -jar \$TACC_PICARD_DIR\/picard\.jar AddOrReplaceReadGroups INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b
# adjust header below...
echo '#!/bin/bash
#SBATCH -J s2b
#SBATCH -n 68
#SBATCH -N 10
#SBATCH -p development
#SBATCH -o cdh.o%j
#SBATCH -e cdh.e%j
#SBATCH -t 1:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
'> heads
# assembling a launching the job:
cat heads s2b >s2b.slurm
sbatch s2b.slurm

rm *sorted*
ls *bt2.bam > bams
cat bams | wc -l  # should be the same number as number of trim files

#-----------
# Calling genotypes using GATK's UnifiedGenotyper

module load gatk
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
echo '#!/bin/bash
#SBATCH -J gt
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -p development
#SBATCH -o gt.o%j
#SBATCH -e gt.e%j
#SBATCH -t 2:00:00
#SBATCH -A [your allocation]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R $GENOME_REF -nt 12 -nct 1 \
--genotype_likelihoods_model SNP \' >unig2
cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unig2
echo '-o primary.vcf ' >> unig2
sbatch unig2

# renaming samples in the vcf file, to get rid of trim-shmim etc
cat primary.vcf | perl -pe 's/\.fq\.trim\.bt2//g' | perl -pe 's/\.trim\.bt2//g' | perl -pe 's/^chrom/chr/' >primary_n.vcf

#----------
# Variant quality score recalibration (VQSR)

# making a tab-delimited table of clone (replicate) sample pairs
nano clonepairs.tab
# paste clone pairs, tab delimited, one pair per line
# Ctl-O , enter, Ctl-X

# extracting "true snps" subset (reproducible across replicates) 
# parameter hetPairs can vary depending on replication scheme (3 is good when you have triplicates)
replicatesMatch.pl vcf=primary_n.vcf replicates=clonepairs.tab hetPairs=3 max.het=0.5 > vqsr.vcf

# determining transition-transversion ratio for true snps (will need it for tranche calibration)
vcftools --vcf vqsr.vcf --TsTv-summary
# Ts/Tv ratio: 1.43  # put your actual number into the next code chunk, --target_titv

# creating recalibration models
module load gatk
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T VariantRecalibrator \
-R $GENOME_REF -input primary_n.vcf -nt 12 \
-resource:repmatch,known=true,training=true,truth=true,prior=30  vqsr.vcf \
-an QD -an MQ -an FS -an HaplotypeScore -mode SNP --maxGaussians 6 \
--target_titv 1.43 -tranche 85.0 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile primary_n.recal -tranchesFile recalibrate.tranches -rscriptFile recalibrate_plots.R 

# examine output and recalibrate*.pdf files - see which tranche to choose the one before TsTv dropoff
# next chunk assumes we are choosing tranche 95

# applying recalibration
module load gatk
export GENOME_REF=where/genome/is/mygenome_ccn.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $GENOME_REF -input primary_n.vcf -nt 12 \
--ts_filter_level 95.0 -mode SNP \
-recalFile primary_n.recal -tranchesFile recalibrate.tranches -o recal.vcf

#---------------
# Applying filters

# identifying poorly genotyped individuals
vcftools --vcf recal.vcf --het
# look at number of sites genotyped per individual (4th column): 
cat out.het 
# see if some samples are much lower in the number of sites than others
# for example, if you want to remove samples showing less than 40000 sites:
cat out.het | awk '$4<40000' | cut -f 1  > underSequenced
cat underSequenced

# applying filter and selecting polymorphic biallelic loci genotyped in 95% or more individuals
# (harsh genotyping rate cutoff is strongly recommended for best quality and to avoid RAD loci affected by null 
# alleles because of mutations in restriction site)
vcftools --vcf recal.vcf --remove underSequenced --remove-filtered-all --max-missing 0.95  --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out filt

# selecting only polymorphic sites (they all are in denovo pipeline!) and sites with no excess heterozygosity
grep -E "#|0/1|0/0.+1/1|1/1.+0/0" filt.recode.vcf >polymorphs.vcf
hetfilter.pl vcf=polymorphs.vcf maxhet=0.5 >best.vcf

#---------------
# Final touches

# genotypic match between pairs of replicates (the most important one is the last one, HetsDiscoveryRate)	
repMatchStats.pl vcf=best.vcf replicates=clonepairs.tab 

# looking at per-individual inbreeding 
# positive - excess false homozygotes (due to poor coverage); negative - false heterozygotes (possibly lumped paralogs)
vcftools --vcf best.vcf --het
cat out.het

# creata a file listing clones and low-site/high-homozygosity individuals to remove
# for example:
echo 'K210
K212
K216
K219
'>clones2remove

# removing clones and badly genotyped ones
vcftools --vcf best.vcf --remove clones2remove --recode --recode-INFO-all --out final

# reinstating original RAD loci (takes about 5 minutes so better launch is as a job or on a idev node
retabvcf.pl vcf=final.recode.vcf tab=cdh_alltags_cc.tab >final.vcf

# thinning for Fst / PCA / ADMIXTURE  (choosing one SNP per tag with max allele frequency):
thinner.pl vcf=final.vcf criterion=maxAF >thinMaxaf.vcf

#######################################
#   DONE, move to 2brad_analysis.txt
#######################################


