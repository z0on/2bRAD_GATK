#!/usr/bin/env python
# Kyle Hernandez
# GetHighQualVcfs.py - Takes a multi-sample VCF file and creates N VCF files (N = # samples)
#                      of high quality SNPs

import argparse
import os
import re
import time
import sys
import logging
from scipy import stats

###############################################################################
# Class for VariantCallFile

class VariantCallFile(object):
    """Asserts the VCF file is version VCFv4.1"""
    def __init__(self, handle):
	self.handle  = handle   # Input VCF File
        self.meta,\
        self.header,\
        self.samples = self.get_meta(handle)
        self.cols    = []       # List container for columns

    # Record iterator
    def __iter__(self):
	with open(self.handle, 'rU') as f:
            for line in f:
                if not line.startswith('#'):
		    self.cols = line.rstrip().split('\t')
		    yield self

    # Initializer to get metadata and sample data
    def get_meta(self, handle):
        meta = []
        meta_dict = {}
	with open(handle, 'rU') as f:
            for line in f:
                if line.startswith('##'):
                    meta.append(line)
                elif line.startswith('#CHROM'):
                    return meta, line, line.rstrip().split('\t')[9::]
                else:
                    break

    # Member functions
    def is_variant(self):
        """BOOL: Is the current record a variant?"""
        if self.ALT() != '.':
            return True
        return False

    def write_row(self, o):
        """Writes out a row of a VCF file in the correct format"""
        o.write('\t'.join(self.cols) + '\n')

    def write_sample_row(self, o, sample_call):
        """
         Writes out a row for a single-indivual VCF file from a 
         multi-individuals VCF file.
        """
        curr_dat = self.cols[:9]
        o.write('\t'.join(curr_dat) + '\t' + sample_call + '\n')
     
    ##########################################################
    # Getters for each element of the VCF file
    def CHROM(self): 
        try: return int(self.cols[0])
	except: return self.cols[0]

    def POS(self): return int(self.cols[1])

    def ID(self): return self.cols[2]
    
    def REF(self): return self.cols[3]

    def ALT(self): return self.cols[4]

    def QUAL(self): return float(self.cols[5])

    def FILTER(self): return self.cols[6]

    def INFO(self): return self.cols[7]

    def FORMAT(self): return self.cols[8]

    def CALLS(self): return self.cols[9::] # List of genotype calls of size Nsamp 

###############################################################################
# Main application
#

def main():
    """
    Main function wrapper.
    """
    logger.info("Initializing VCF file...")
    # Initialize the VCF object
    vcf_init = VariantCallFile(args.infile)
    
    # Estimate cutoff
    logger.info("Estimating the " + str(args.percentile) + "th percentile cutoff...")
    cutoff = extract_quals(vcf_init, args.percentile)
    logger.info("The cutoff is " + str(cutoff))

    # Filter VCF
    logger.info("Filtering VCF file...")
    filter_file(vcf_init, cutoff, args)

def filter_file(vcf_records, cutoff, args):
    """
    Writes variants with GQ > args.GQ and MAPQ > percentile cutoff
    """
    flag       = 0
    parent_dir = os.path.abspath(args.outdir) + os.sep
    fil_list   = []
    header     = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    for record in vcf_records:
        # First open all the files for single/multi-sampled VCF.
        # Then write the meta-data and column header with correct sample name.
        if flag == 0:
            fil_list = [open(parent_dir + i + "_HQ.vcf", 'wb') for i in record.samples]
            logger.info("Writing high quality SNPs to " + str(len(fil_list)) + " individual VCFs...")
            [j.write(''.join(record.meta)) for j in fil_list]
            [j.write(header + '\t' + record.samples[n] + '\n') for n,j in enumerate(fil_list)]
            flag = 1
        # Next, we filter out positions with low MAPQ
        if record.is_variant() and record.QUAL() > cutoff:
            # Here we need to then only print out snps to samples where
            # GQ > the cutoff
            # But we only want to print out samples where GT != 0/0 which means it's the reference
            if args.ploidy == 1:
                for n,i in enumerate(record.CALLS()):
                    try:
                        if ':' in i and int(i.split(':')[3]) > args.GQ and i.split(':')[0] != '0':
                            record.write_sample_row(fil_list[n], i)
                    except IndexError:
                        logger.warn(
                          "Weird genotype record '{0}' at CHROM '{1}' POS '{2}' SAMPLE '{3}'".format(
                          i, record.CHROM(), record.POS(), record.samples[n]))
                        pass 
            else:
                for n,i in enumerate(record.CALLS()):
                    try:
                        if ':' in i and int(i.split(':')[3]) > args.GQ and i.split(':')[0] != '0/0':
                            record.write_sample_row(fil_list[n], i)
                    except IndexError:
                        logger.warn(
                          "Weird genotype record '{0}' at CHROM '{1}' POS '{2}' SAMPLE '{3}'".format(
                          i, record.CHROM(), record.POS(), record.samples[n]))
                        pass 
    [j.close() for j in fil_list]

def extract_quals(vcf_records, p):
    """
    Adapted from J. Malcom
    Extract ALTQs for variant loci in vcf record; return percentile score
    """
    AQ = [rec.QUAL() for rec in vcf_records if rec.ALT() != '.']
    return stats.scoreatpercentile(AQ, p)

if __name__ == '__main__':
    start = time.time()

    # Initialize logger
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Main parser
    parser = argparse.ArgumentParser(prog='GetHighQualVcfs.py',\
	description = 'Split multi-sample VCFs into single sample VCFs of high quality SNPs.')

    # Command line args
    parser.add_argument('-i', '--infile', required=True, type=str,\
	help='Multi-sample VCF file')
    parser.add_argument('-o', '--outdir', required=True, type=str,\
	help='Directory to output HQ VCF files.')
    parser.add_argument('--ploidy', type=int, default=2,\
	help='1 for haploid; 2 for diploid')
    parser.add_argument('--GQ', default=90, type=int,\
	help='Filters out variants with GQ < this limit.')
    parser.add_argument('--percentile', default=90, type=int,\
	help='Reduces to variants with ALTQ > this percentile.')

    # Initialize the parser
    args = parser.parse_args()
   
    # Run script
    logger.info('-'*80)
    logger.info('Kyle Hernandez, 2013, kmhernan84@gmail.com')
    logger.info('GetHighQualVcfs - Get high quality SNPs from a multi-sample VCF files for GATK') 
    logger.info('-'*80)
    main()
     
    logger.info("Finished; Took: " + str(time.time() - start) + " seconds.")
