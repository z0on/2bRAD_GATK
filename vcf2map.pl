#!/usr/bin/env perl

my $usage="

vcf2map.pl : replaces scaffold coordinates in the vcf file with their approximate locations 
in the linkage map. Scaffolds not pinned by the map are discarded.
The input VCF file is assumed to be sorted by scaffold.

Arguments:

vcf=[filename]
map=[filename] : table where third column is linkage group, fourth column is map position,
                 fifth column is the scaffold ID, and sixth column is scaffold position.
     mbcm=0.33 : megabases per centimorgan 

Outputs two files mapAnchored_[input vcf filename] :
.vcf: recoded sorted vcf
.tab: table of old to new coordinates:
      scaffold	scaff.coordinate	chromosome	chrom.coordinate

Mikhail Matz, matz\@utexas.edu

";

my $vcf="";
my $mapfile="";
my $mbcm=0.33;
if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage;}
if ("@ARGV"=~/map=(\S+)/) { $mapfile=$1; } else { die $usage;}
if ("@ARGV"=~/mbcm=(\S+)/) { $mbcm=$1; }

open MAP, $mapfile or die "cannot open map table $map\n";
my %map={};
my %lg={};
my %coord={};

while (<MAP>){
	chop;
	my @line=split("\t",$_);
	$lg{$line[4]}=$line[2];
	if (" @{$map{$line[4]}} "!~/$line[3]/) {
		push @{$lg{$line[4]}},$line[2];
		push @{$map{$line[4]}}, $line[3] ;
		push @{$coord{$line[4]}}, $line[5];
	}
}
close MAP;

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my $s0="";
my @c=();
my @g=();
my %lgAdd={};
while (<VCF>) {
	(my $s,my $co,my @rest)=split(/\t/,$_);	
	next unless @{$map{$s}}>0;
	my $nearest=1000000000000;
	my $spot=0;
	for (my $i=0; my $c=${$coord{$s}}[$i];$i++){
		my $dist=$co-$c;
		if (abs($dist)<$nearest) {
			$spot=$i;
			$nearest=abs($dist);
		}
	}
	my $loc=$mbcm*(10**6)*${$map{$s}}[$spot]+($co-${$coord{$s}}[$spot]);
	if ($loc<0 and abs($loc)>$lgAdd{${$lg{$s}}[$spot]}) { 
		$lgAdd{${$lg{$s}}[$spot]}=abs($loc)+1; 
#warn "lg:${$lg{$s}}[$spot] add:$lgAdd{${$lg{$s}}[$spot]}\n";
	}
}

close VCF;
open VCF, $vcf;
my $outvcf="mapAnchored_".$vcf;
open OUTV, ">$outvcf" or die "cannot create output vcf file $outvcf\n";

my %data={};
my %tab={};
while (<VCF>) {
	if ($_=~/^#/) { print {OUTV} $_ and next;}
	(my $s,my $co,my @rest)=split(/\t/,$_);	
	next unless @{$map{$s}}>0;
	my $nearest=1000000000000;
	my $spot=0;
	for (my $i=0; my $c=${$coord{$s}}[$i];$i++){
		my $dist=$co-$c;
		if (abs($dist)<$nearest) {
			$spot=$i;
			$nearest=abs($dist);
		}
	}
	my $loc=$mbcm*(10**6)*${$map{$s}}[$spot]+($co-${$coord{$s}}[$spot])+$lgAdd{${$lg{$s}}[$spot]};
#print "$s $co|@{$coord{$s}}|nearest:$spot($nearest) lg:${$lg{$s}}[$spot] map:${$map{$s}}[$spot] loc:$loc\n";
	my $ind=(10**9)*${$lg{$s}}[$spot]+$loc;
	$data{$ind}=join("\t",@rest);
	$data{$ind}="chr${$lg{$s}}[$spot]\t$loc\t".$data{$ind};
	$tab{$ind}="$s\t$co\tchr${$lg{$s}}[$spot]\t$loc\n";
}


my $outtab=$outvcf;
$outtab=~s/\.vcf/\.tab/;
open OUTT, ">$outtab" or die "cannot create output tab file $outtab\n";

my @indx=sort {$a <=> $b } keys %data;
my $count=0;
foreach my $i (@indx) {  
	next if ($i=~/HASH/);
	print {OUTV} $data{$i};
	print {OUTT} $tab{$i};
	$count++;
}
warn "\n$count variants anchored\n\n";
	