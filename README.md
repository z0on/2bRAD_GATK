Whole genome reference-based genotyping with 2bRAD
------------------------------------------

Mikhail Matz, matz@utexas.edu

2bRAD has been described in Wang et al 2012 
http://www.nature.com/nmeth/journal/v9/n8/abs/nmeth.2023.html 

The advanced features of 2bRAD include removal of PCR duplicates, as well as advanced filtering and assessment of overall genotyping quality based on genotyping replicates
This collection of scripts is for using 2bRAD data with GATK pipeline

This repository contains the lab protocol for sample preparation, as well as scripts and walkthrough (2bRAD_gatk_README.txt) for:
- trimming and quality filtering;
- removing PCR duplicates;
- mapping to reference genome;
- recalibrating base qualities and realigning around indels;
- varint calling using UnifiedGenotyper
- recalibrating variant quality scores based on genotyping replicates;
- smart-thinning and final filtering;
- quality assessment based on replicates.
