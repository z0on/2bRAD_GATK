GATK pipeline for 2bRAD data, based on bowtie2 mapping to a reference genome
Last updated November 20, 2017, Mikhail Matz matz@utexas.edu (includes ANGSD)

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
PGDspider2: http://www.cmpg.unibe.ch/software/PGDSpider/

see installationInstructions_3dPartySoftware.txt for details on installing ANGSD and some other things

NOTE: this walkthrough has been written for lonestar cluster of the Texas
Advanced Computer Center, which runs 48 tasks per node and uses SLURM batch 
processing system. To adopt the walkthough to your cluster you 
might need to edit the launcher_creator.py script to make it
compatible with your cluster. 

==============================================
IMPORTANT: DO NOT SEQUENCE 2BRAD LIBRARIES ALONE ON A HISEQ4000 LANE! INVARIANT BASES
(ADAPTOR, RESTRICTION SITE) WILL NOT BE READ WELL. MIX IN 20% OF PHIX SAMPLES TO 
DIVERSIFY; OR MIX YOUR 2BRAD LIBRARIES WITH SOMEONE ELSE'S NOT-2BRAD LIBRARIES 
(JUST MAKE SURE THE BARCODES DO NOT OVERLAP) 

If you did and you have used BcgI enzyme, use 2bRAD_trim_launch_dedup2.pl for trimming
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
launcher_creator.py -j trims -n trims -l trimjob -t 0:30:00 -a mega2014 -e matz@utexas.edu
sbatch trimjob

# NB: use this command instead of the one above if you have 2bRAD libraries without 
# degenerate tag but with 4-base in-line barcode:
2bRAD_trim_launch.pl fastq barcode2=4 > trims
# And if you have the original-style libraries without degenerate tag and without in-line barcode, use this:
2bRAD_trim_launch.pl fastq > trims
# AND IF you sequenced your double-barcoded, deduplicatable libraries on HiSeq 4000 alone
# (resulting in poor quality at restriction site and adaptor bases) and you have used BcgI enzyme, use this:
2bRAD_trim_launch_dedup2.pl fastq > trims

# quality filtering using fastx_toolkit
module load fastx_toolkit
ls *.tr0 | perl -pe 's/^(\S+)\.tr0$/cat $1\.tr0 \| fastq_quality_filter -q 20 -p 90 >$1\.trim/' >filt0

# NOTE: run the next line ONLY if your qualities are 33-based 
# (if you don't know just try to see if it works eventually, if you get errors from fastx_toolkit, try the other one):
	cat filt0 | perl -pe 's/filter /filter -Q33 /' > filt
#if you did NOT run the line above, run this one (after removing # symbol):
#	mv filt0 filt

launcher_creator.py -j filt -n filt -l filtjob -t 0:30:00 -a mega2014 -e matz@utexas.edu
sbatch filtjob

#----------------------------------
# setting up the genome reference

# concatenating genome contigs into small number of "pseudo-chromosomes" (to help memory usage in GATK)
concatFasta.pl fasta=mygenome.fasta

# normalize fasta records to ensure all lines are of the same length, using Picard
module load picard
java -Xmx1g -jar $TACC_PICARD_DIR/NormalizeFasta.jar INPUT=mygenome_cc.fasta OUTPUT=mygenome_ccn.fasta

# Now, edit these lines THROUGHOUT THE TEXT to fit the location and name of your genome reference 
export GENOME_FASTA=mygenome_ccn.fasta
export GENOME_DICT=mygenome_ccn.dict 

bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA
java -jar $TACC_PICARD_DIR/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

#-------------------------------------
# mapping to the reference genome and creating BAM files
export GENOME_FASTA=mygenome_ccn.fasta

# creating list of commands out of file list using 'perl -pe'
ls *trim | perl -pe 's/(\S+)/bowtie2 --no-unal -x \$GENOME_FASTA -U $1 -S $1\.sam/' > bt2
launcher_creator.py -j bt2 -n maps -t 2:00:00 -a mega2014 -e matz@utexas.edu
sbatch maps.slurm

# mapping efficiency?
tail -100 maps.e*

# what does this line do?
ls *.sam > sams

# next stage is compressing, sorting and indexing the SAM files (so they become BAM files)
export GENOME_FASTA=mygenome_ccn.fasta
cat sams | perl -pe 's/(\S+)\.sam/samtools import \$GENOME_FASTA $1\.sam $1\.unsorted\.bam && samtools sort -o $1\.sorted\.bam $1\.unsorted\.bam && java -Xmx5g -jar \$TACC_PICARD_DIR\/picard\.jar AddOrReplaceReadGroups INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b
launcher_creator.py -j s2b -n s2b -t 0:05:00 -a mega2014 -e matz@utexas.edu
sbatch s2b.slurm

rm *sorted*
ls *bam | wc -l  # should be the same number as number of trim files

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis - save them for your own data

#==========================
# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).

# open interactive session
idev

# Generating genotype likelihoods from highly confident (non-sequencing-error) SNPs
# ANGSD parameters to use:

# ALWAYS record your filter settings and explore different combinations to confirm that results are robust. 
# Our filters will be:
# -minMapQ 30 only highly unique mappings
# -minQ 30 only highly confident base calls
# -minInd 50 the site must be genotyped in at least 50 individuals (set this to at least 80% of your total number of your individuals) 
# -snp_pval 1e-5 high confidence that the SNP is not just sequencing error 
# -minMaf 0.05 only common SNPs, with allele frequency 0.05 or more.
# Note: the last two filters are very efficient against sequencing errors but introduce bias against true rare alleles. It is OK (and even desirable) - UNLESS we want to do AFS analysis. We will generate data for AFS analysis in the next part.
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 50 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5"

# things we want ANGSD to do: 
# -GL 1 : samtools likelihood model
# -doGlf 2 : output beagle format (for admixture)
# -doGeno 32 : binary genotype likelihoods format (for ngsCovar => PCA)
# -doMajorMinor 1 : infer major and minor alleles from data (not from reference)
DOS="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2"
angsd -b bams -GL 1 $FILTERS $DOS -P 40 -out commonSNPs

# how many SNPs?
NSITES=`zcat commonSNPs.beagle.gz | wc -l`
echo $NSITES

# covariance matrix for PCA (if coverage is approximately equal across samples):
gunzip commonSNPs.geno.gz
ngsCovar -probfile commonSNPs.geno -outfile commonSNPs.covar -nind 61 -nsites $NSITES -call 0 -norm 0 

# if coverage is very low and/or unequal, use commonSNPs.covMat and commonSNPs.ibsMat from angsd run for PCoA and PCA. In fact, using these results would be conservative in any case.

# ADMIXTURE for K from 2 to 5
for K in `seq 2 5` ; 
do 
NGSadmix -likes commonSNPs.beagle.gz -K $K -P 10 -o commonSNPs_k${K};
done

# scp *Mat, *covar, *qopt and bams files to laptop, use angsd_ibs_pca.R and admixturePlotting_v4.R to plot PCA and ADMIXTURE

#==========================
# ANDSD => SFS for demographic analysis

# first, manually edit file bams in nano to remove all genotyping replicates (if you had them)
# enter interactive session
idev

# creating list of SNP sites for SFS, for all pops:
# running ANGSD without snp_pval and  minMaf filter 
# using stronger -minInd filter (set it to 90% of total number of individuals)

FILTERS="-minMapQ 30 -minQ 30 -minInd 50 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5"
 DOS="-doMajorMinor 1 -doMaf 1"
angsd -b bams_noclones -GL 1 $FILTERS $DOS -P 48 -out sfsSites

# extracting and sorting list of sites to make the SFS out of for all pops
zcat sfsSites.mafs.gz | perl -pe 's/chr//'| cut -f 1,2 | tail -n +2 | sort -k 1,1n -k 2,2n | perl -pe 's/\s/:/' | perl -pe 's/^(.)/chr$1/' > sites2do

# exit interactive session
exit

# how many sites?
cat sites2do | wc -l

# making lists of bam files for each pop (assuming name of population - in this case, KIE, RAU and ESP - was a part of read file name)
grep KIE bams_noclones >kie
grep RAU bams_noclones >rau
grep ESP bams_noclones >esp

# the next job takes quite a bit of time, ~1h
echo 'angsd -b kie -rf sites2do -GL 1 -doSaf 1 -anc mygenome_ccn.fasta -P 1 -out kie
angsd -b esp -rf sites2do -GL 1 -doSaf 1 -anc mygenome_ccn.fasta -P 1 -out esp
angsd -b rau -rf sites2do -GL 1 -doSaf 1 -anc mygenome_ccn.fasta -P 1 -out rau' >sfs0
launcher_creator.py -j sfs0 -n sfs0 -t 2:00:00 -a mega2014 -e matz@utexas.edu
sbatch sfs0.slurm

# exporting probabilities of allele counts
realSFS print kie.saf.idx > kie.sfs
realSFS print esp.saf.idx > esp.sfs
realSFS print rau.saf.idx > rau.sfs

# making dadi-SNP counts table based on genotype probabilities
ls *.sfs >sfss
Rscript ~/bin/sfs2dadi.R infiles=sfss  prefix=commonSNPs

# exporting folded 1d-sfs out of dadi-SNP format
dadi2sfs_fold.py commonSNPs_dadi.data kie 40
dadi2sfs_fold.py commonSNPs_dadi.data esp 40
dadi2sfs_fold.py commonSNPs_dadi.data rau 42

#-----------------
# 2d AFS analysis using moments

# install moments (python package, see the beginning of this file)

# get Misha's moments scripts collection
cd
git clone https://github.com/z0on/AFS-analysis-with-moments.git
chmod +x AFS-analysis-with-moments/*
mv AFS-analysis-with-moments/* ~/bin
rm -rf AFS-analysis-with-moments 
cds
cd RAD

# NB: must mask singletons and private alleles during SFS analysis! 

# interactive session
idev

# print 2d SFS:
2dAFS_fold.py commonSNPs_dadi.data kie rau 40 42
2dAFS_fold.py commonSNPs_dadi.data esp rau 40 42
2dAFS_fold.py commonSNPs_dadi.data kie esp 40 40
2dAFS.py commonSNPs_dadi.data kie rau 40 42
2dAFS.py commonSNPs_dadi.data esp rau 40 42
2dAFS.py commonSNPs_dadi.data kie esp 40 40

#scp *.pdf files to laptop, look at them

# S2M model for kie and rau, masking singletons (nu1, nu2, T, m12, m21 = params):
# paste this same command 20 times in the text file in nano, create and launch parallel job from it using launcher_creator.py (as before). 
S2M_fold_masked.py commonSNPs_dadi.data kie rau 40 42 1 1 1 1 1

# IM2 model for kie and rau (nu1_0,nu2_0,nu1,nu2,T,m12,m21 = params):
IM2_fold_masked.py commonSNPs_dadi.data kie rau 40 42 4 7 4 7 13 2.5 0.7

#-----------------

# stairway plot: historical population sizes for ESP
# (must install stairwayPlot first, and replace $WORK/stairway_plot_v2beta/stairway_plot_es below with your own path to stairway_plot files

(ignoring singletons)

realSFS esp.saf.idx >esp_sfs
cat esp_sfs
3759.234804 5082.773410 4901.517283 2725.871576 1952.031534 1325.980506 1176.750821 941.780456 567.845373 817.837276 401.409512 984.552254 16.381838 644.826881 61.597748 565.959387 2.839908 442.083707 94.392939 113.486667 218.516297 82.584333 87.582359 11.942078 0.002080 69.575571 4.658302 0.000226 0.000962 43.606802 0.079898 0.001061 10.795316 0.000000 24.617357 4.025453 47.096118 37.027363 59.100985 69.955557 50.678000 

cat esp_sfs_folded
41 
6537 6078 4625 3291 2351 1800 1347 1094 949 810 734 641 573 491 436 384 418 364 343 336 189 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

# copy non-zero entries except the first one 
nano esp.blueprint
#example blueprint file
#input setting
popid: ESP # id of the population (no white space)
nseq: 40 # number of sequences
L: 700000 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whethr the SFS is folded (true or false)
SFS: 6078 4625 3291 2351 1800 1347 1094 949 810 734 641 573 491 436 384 418 364 343 336 189
#smallest_size_of_SFS_bin_used_for_estimation: 2 # default is 1; to ignore singletons, change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: # default is n-2; to ignore singletons, change this number to nseq-2
pct_training: 0.67 # percentage of sites for training
nrand: 9 18 28 38 # number of random break points for each try (separated by white space)roughly (nseq-2)/4, (nseq-2)/2, (nseq-2)*3/4, nseq-2
project_dir: ESP # project directory
stairway_plot_dir: $WORK/stairway_plot_v2beta/stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#output setting
mu: 2e-8 # assumed mutation rate per site per generation
year_per_generation: 0.25 # assumed generation time (in years)
#plot setting
plot_title: esp # title of the plot
xrange: 0.1,1000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size


java -cp $WORK/stairway_plot_v2beta/stairway_plot_es Stairbuilder esp.blueprint 
grep swarmops esp.blueprint.sh > addTheta
launcher_creator.py -j addTheta -n at -e matz@utexas.edu -a tagmap -t 0:15:00 -q normal

nano at.slurm
# add:
#SBATCH -N 18

sbatch at.slurm

grep -v swarmops esp.blueprint.sh >movesPlots
bash movesPlots

scp ESP/*.summary* to laptop

#==========================

#        G  A  T  K

# ("hard-call" genotyping, use only for high-coverage data, >10x after deduplication)

export GENOME_REF=mygenome_ccn.fasta
ls *.bam > bams

# writing command script
echo '#!/bin/bash
#SBATCH -J gt
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p development
#SBATCH -o gt.o%j
#SBATCH -e gt.e%j
#SBATCH -t 2:00:00
#SBATCH -A mega2014
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matz@utexas.edu
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R $GENOME_REF -nt 40 -nct 1 \
--genotype_likelihoods_model SNP \' >unig2
cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unig2
echo '-o primary.vcf ' >> unig2
sbatch unig2

#----------
# Variant quality score recalibration (VQSR)

# making a tab-delimited table of clone (replicate) sample pairs
nano clonepairs.tab
# paste clone pairs, tab delimited, one pair per line, like this:
IB14ESP07A.trim	IB14ESP07B.trim
IB14ESP09A.trim	IB14ESP09B.trim
IB14ESP21A.trim	IB14ESP21B.trim
IB14ESP30A.trim	IB14ESP30B.trim
IB14KIE04a.trim	IB14KIE04b.trim
IB14KIE07a.trim	IB14KIE07b.trim
IB14KIE08a.trim	IB14KIE08b.trim
IB14KIE13a.trim	IB14KIE13b.trim
IB14RAU05a.trim	IB14RAU05b.trim
# Ctl-O , enter, Ctl-X

# extracting "true snps" subset (reproducible across replicates) 
# parameter hetPairs can vary depending on replication scheme (3 is good when you have triplicates)
replicatesMatch.pl vcf=primary.vcf replicates=clonepairs.tab hetPairs=2 max.het=0.5 > vqsr.vcf

# determining transition-transversion ratio for true snps (will need it for tranche calibration)
vcftools --vcf vqsr.vcf --TsTv-summary
# Ts/Tv ratio: 1.244  # put your actual number into the next code chunk, --target_titv

# creating recalibration models
export GENOME_REF=mygenome_ccn.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T VariantRecalibrator \
-R $GENOME_REF -input primary_n.vcf -nt 12 \
-resource:repmatch,known=true,training=true,truth=true,prior=30  vqsr.vcf \
-an QD -an MQ -an FS -mode SNP --maxGaussians 6 \
--target_titv 1.244 -tranche 85.0 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile primary.recal -tranchesFile recalibrate.tranches -rscriptFile recalibrate_plots.R 

# examine output and recalibrate*.pdf files - see which tranche to choose the one before TsTv dropoff
# next chunk assumes we are choosing tranche 95

# applying recalibration (95% tranche)
export GENOME_REF=mygenome_ccn.fasta
java -jar $TACC_GATK_DIR/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $GENOME_REF -input primary_n.vcf -nt 12 \
--ts_filter_level 95.0 -mode SNP \
-recalFile primary.recal -tranchesFile recalibrate.tranches -o recal.vcf

#---------------
# Applying filters

# Tjarno peeps: SKIP THIS BIT
# identifying poorly genotyped individuals
vcftools --vcf recal.vcf --het
# look at number of sites genotyped per individual (4th column): 
cat out.het 
# see if some samples are much lower in the number of sites than others
# for example, if you want to remove samples showing less than 40000 sites:
cat out.het | awk '$4<40000' | cut -f 1  > underSequenced
cat underSequenced

# applying filter and selecting polymorphic biallelic loci genotyped in 90% or more individuals
# (harsh genotyping rate cutoff is strongly recommended for best quality and to avoid RAD loci affected by null 
# alleles because of mutations in restriction site)
vcftools --vcf recal.vcf --remove-filtered-all --max-missing 0.9  --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out filt

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

# create a file listing clones (and low-site/high-homozygosity individuals, if any) to remove
cat clonepairs.tab | cut -f 2 >clones2remove

# removing clones and badly genotyped ones
vcftools --vcf best.vcf --remove clones2remove --recode --recode-INFO-all --out final

# thinning for Fst / PCA / ADMIXTURE  (choosing one SNP per tag with max allele frequency):
thinner.pl vcf=final.recode.vcf criterion=maxAF >thinMaxaf.vcf

#----------
# ADMIXTURE using hard-call data (vcf)

# creating a dataset with fake "chr" chromosome designations
cat thinMaxaf.vcf | perl -pe 's/tag(\d)(\d+)\t(\d+)/chr$1\t$3$2/'>thinMaxChrom.vcf

# reformatting VCF into plink binary BED format
plink --vcf thinMaxChrom.vcf --make-bed --out rads

# ADMIXTURE with cross-validation to select K 
# (bash script to run admixture with several different K's)
for K in 1 2 3 4 5; do admixture --cv rads.bed $K | tee log${K}.out; done

# minimal cross-validation error = optimal K
grep -h CV log*.out

grep '#CHROM' thinMaxChrom.vcf | perl -pe 's/\t/\n/g' | tail -n +10 | perl -pe 's/^IB14(...)(.+)/IB14$1$2\t$1/' > inds2pops

# scp the *.Q and inds2pops to laptop, plot it in R:
tbl=read.table("rads.3.Q")
barplot(t(as.matrix(tbl)), col=rainbow(5),xlab="Individual #", ylab="Ancestry", border=NA)

# or, more fancy, use admixturePlotting_V4.R

