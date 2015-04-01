#!/usr/bin/perl

my $usage="

filter_stats_gatk.vcf : 

Calculates quantiles for three filtering criteria using VCF file:
- fisher-strand (FS)
- mapping quality (MQ)
- depth (DP)

For MQ, low values are bad and high values are good;
for FS, high values are bad;
for DP, both high and low values are bad

Arguments:

vcf=[file name]  vcf file to analyze

Output: 
Tables of quantiles (printed to STDOUT).

Example:
filter_stats_gatk.pl vcf=round2.vcf

NOTE: only use this to select your filtering criteria if you don't have 
genotyping replicates.

";

my $vcf;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @dp=();
my @fs=();
my @mq=();
my @dline=();
my @iline=();
my $geno;

while (<VCF>) {
	next if ($_=~/^#/);
	@dline=split("\t",$_);
	$geno=()=$_=~/\d\/\d/g;
	next if ($geno==0);
	@iline=split(";",$dline[7]);
	foreach my $i (@iline){
		if ($i=~/MQ=(\d+)/) { push @mq, $1;}
		elsif ($i=~/FS=(\d+)/) { push @fs, $1;}
		elsif ($i=~/DP=(\d+)/) { push @dp, sprintf("%.0f",$1/$geno);}
	}
}
close VCF;

@mq=sort {$a <=> $b} @mq;
@fs=sort {$b <=> $a} @fs;
@dp=sort {$a <=> $b} @dp;

my $total=$#dp+1;

my %mqq;
my %fsq;
my %dpq;
my %dpq2;

$mqq{1}=$mq[sprintf("%.0f",$total*0.01)];
$fsq{1}=$fs[sprintf("%.0f",$total*0.01)];
$mqq{5}=$mq[sprintf("%.0f",$total*0.05)];
$fsq{5}=$fs[sprintf("%.0f",$total*0.05)];
$mqq{10}=$mq[sprintf("%.0f",$total*0.1)];
$fsq{10}=$fs[sprintf("%.0f",$total*0.1)];
$mqq{15}=$mq[sprintf("%.0f",$total*0.15)];
$fsq{15}=$fs[sprintf("%.0f",$total*0.15)];
$mqq{20}=$mq[sprintf("%.0f",$total*0.2)];
$fsq{20}=$fs[sprintf("%.0f",$total*0.2)];
$mqq{30}=$mq[sprintf("%.0f",$total*0.3)];
$fsq{30}=$fs[sprintf("%.0f",$total*0.3)];
$mqq{40}=$mq[sprintf("%.0f",$total*0.4)];
$fsq{40}=$fs[sprintf("%.0f",$total*0.4)];
$mqq{50}=$mq[sprintf("%.0f",$total*0.5)];
$fsq{50}=$fs[sprintf("%.0f",$total*0.5)];
$mqq{60}=$mq[sprintf("%.0f",$total*0.6)];
$fsq{60}=$fs[sprintf("%.0f",$total*0.6)];
$mqq{70}=$mq[sprintf("%.0f",$total*0.7)];
$fsq{70}=$fs[sprintf("%.0f",$total*0.7)];
$mqq{80}=$mq[sprintf("%.0f",$total*0.8)];
$fsq{80}=$fs[sprintf("%.0f",$total*0.8)];
$mqq{90}=$mq[sprintf("%.0f",$total*0.9)];
$fsq{90}=$fs[sprintf("%.0f",$total*0.9)];


$dpq{1}=$dp[sprintf("%.0f",$total*0.005)];
$dpq2{1}=$dp[sprintf("%.0f",$total*0.995)];
$dpq{5}=$dp[sprintf("%.0f",$total*0.025)];
$dpq2{5}=$dp[sprintf("%.0f",$total*0.975)];
$dpq{10}=$dp[sprintf("%.0f",$total*0.05)];
$dpq2{10}=$dp[sprintf("%.0f",$total*0.95)];
$dpq{15}=$dp[sprintf("%.0f",$total*0.075)];
$dpq2{15}=$dp[sprintf("%.0f",$total*0.925)];
$dpq{20}=$dp[sprintf("%.0f",$total*0.1)];
$dpq2{20}=$dp[sprintf("%.0f",$total*0.9)];
$dpq{30}=$dp[sprintf("%.0f",$total*0.15)];
$dpq2{30}=$dp[sprintf("%.0f",$total*0.85)];
$dpq{40}=$dp[sprintf("%.0f",$total*0.2)];
$dpq2{40}=$dp[sprintf("%.0f",$total*0.8)];
$dpq{50}=$dp[sprintf("%.0f",$total*0.25)];
$dpq2{50}=$dp[sprintf("%.0f",$total*0.75)];
$dpq{60}=$dp[sprintf("%.0f",$total*0.3)];
$dpq2{60}=$dp[sprintf("%.0f",$total*0.7)];
$dpq{70}=$dp[sprintf("%.0f",$total*0.35)];
$dpq2{70}=$dp[sprintf("%.0f",$total*0.65)];
$dpq{80}=$dp[sprintf("%.0f",$total*0.4)];
$dpq2{80}=$dp[sprintf("%.0f",$total*0.6)];
$dpq{90}=$dp[sprintf("%.0f",$total*0.45)];
$dpq2{90}=$dp[sprintf("%.0f",$total*0.55)];

print "\nMQ quantiles (low numbers are bad):\n";
foreach my $q (sort {$a <=> $b} keys %mqq){
	print 100-$q,"%:\t>$mqq{$q}\n";
}

print  "\nFS quantiles (high numbers are bad):\n";
foreach my $q (sort {$a <=> $b} keys %fsq){
	print 100-$q,"%:\t<$fsq{$q}\n";
}

print  "\nmeanDP quantiles (extreme numbers are bad):\n";
foreach my $q (sort {$a <=> $b} keys %dpq){
	print  100-$q,"%:\t$dpq{$q} - $dpq2{$q}\n";
}


