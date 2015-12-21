#!/usr/bin/perl

my $usage="
vcf2fasta_ref.pl

Takes a vcf file and the mapping reference fasta file,
extracts a stretch of DNA sequence surrounding each SNP.

Arguments:

       vcf=[filename] : vcf file name
reference=[file name] : reference file name
     spread=[integer] : number of bases flanking the SNP on each side to extract. 
                        Default is 50 (results in 101 b fragments)

Prints to STDOUT

Mikhail Matz, matz\@utexas.edu

";

use warnings;

my $vcf;
my $ref;
my $spread=50;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }
if ("@ARGV"=~/reference=(\S+)/) { $ref=$1;}
else { die $usage; }
if ("@ARGV"=~/spread=(\d+)/) { $spread=$1;}

# reading reference

warn "\nreading genome....\n";

my %genome={};
open REF, "$ref";
my $seq;
my $head;

while (<REF>) {
	chop;
	if ($_=~/^>(\S+)/){
		if ($seq){
			$genome{$head}=$seq;
		}
		$head=$1;
		$seq="";
	}
	else {
		$seq.=$_;
	}
}
if ($seq){
	$genome{$head}=$seq;
}

close REF;

#my @heads=keys(%genome);

open VCF, "$vcf";

while (<VCF>) {
	next if ($_=~/^#/);
	chomp;
	(my $chr, my $pos, my @rest)=split('/\t/',$_);
	if ($genome[$chr] && length($genome[$chr])>$pos+$spread) {
		print ">$chr|$pos|$rest[1]$rest[2]\n", substr($genome[$chr],$pos-$spread,2*$spread+1),"\n";
	}
	else { warn "cannot extract $chr $pos\n"; }
}
close VCF;
	
	
	

