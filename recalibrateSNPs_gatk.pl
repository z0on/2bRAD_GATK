#!/usr/bin/perl

my $usage="

recalibrateSNPs_gatk.pl  :

Non-parametric variant quality recalibration based on INFO fields 
in the set of variants that are reproducibly genotyped among replicates 
(output of replicatesMatch.pl)

This script is for reference-based pipeline. Use recalibrateSNPs.pl for
de novo pipeline.

Three parameters (INFO field) are assessed for each SNP: 
- fisher-strand (FS)
- mapping quality (MQ)
- quality by depth (QD) 
- depth (DP)

The script computes the product of quantiles for different INFO fields, determines 
the quantiles of the result within the 'true' set (JOINT quantiles),
and then computes the new quality scores in the main vcf file.

These scores are supposed to correspond to the probability (x100) that the 
SNPs comes from the same distribution as the 'true' SNPs. 

Output:

- Recalibrated VCF (printed to STDOUT) with QUAL field replaced by the new score 
  minus 0.1 (for easy filtering)
  
- a table (printed to STDERR) showing the \"gain\" at each recalibrated quality 
  score setting, which is excess variants removed when filtering at this quality score.
  For example, if 40% of all variants are removed at the quality score 15, the gain 
  is 40 - 15 = 35% (i.e., in addition to removing 15% of the 'true' variants, 
  additional 35% of the dataset, likely corresponding to wrong variants, is removed). 
  The optimal filtering score is the one giving maximum gain. Try different combinations 
  of the four possible filters to find the one maximizing the gain.

Arguments:

vcf=[file name]  : vcf file to be recalibrated

true=[file name] : vcf file that is subset of the above, with SNPs that are considered 
                   true because they show matching and polymorphic genotypes in replicates 
                   (replicatesMatch.pl polyonly=1). 
                   
           -nodp : do not use DP 
           -nofs : do not use FS
           -nomq : do not use MQ
           -noqd : do not use QD

Example: 
recalibrateSNPs_gatk.pl vcf=round2.names.vcf true=vqsr.vcf > recal.vcf

Mikhail Matz, matz\@utexas.edu September 30, 2014
 
";

my $vcf;
my $true;
my $nodp=0;
my $nofs=0;
my $nomq=0;
my $noqd=0;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }
if ("@ARGV"=~/true=(\S+)/) { $true=$1;}
else { die $usage; }
if ("@ARGV"=~/-nomq/) { $nomq=1;}
if ("@ARGV"=~/-nofs/) { $nofs=1;}
if ("@ARGV"=~/-nodp/) { $nodp=1;}
if ("@ARGV"=~/-noqd/) { $noqd=1;}
if($nodp+$nofs+$nomq+$noqd==4) { die "no fields left to filter!\n";}

open TR, $true or die "cannot open the true set $true\n";

my @dp=();
my @mq=();
my @fs=();
my @qd=();
my @dline=();
my @iline=();

while (<TR>) {
	next if ($_=~/^#/);
	@dline=split("\t",$_);
	@iline=split(";",$dline[7]);
	foreach my $i (@iline){
		if ($i=~/MQ=([-\d\.]+)/) { push @mq, $1;}
		elsif ($i=~/FS=([-\d\.]+)/) { push @fs, $1;}
		elsif ($i=~/DP=([-\d\.]+)/) { push @dp, $1;}
		elsif ($i=~/QD=([-\d\.]+)/) { push @qd, $1;}
	}
}
close TR;

@mq=sort {$a <=> $b} @mq;
@fs=sort {$b <=> $a} @fs;
@dp=sort {$a <=> $b} @dp;
@qd=sort {$b <=> $a} @qd;

my $total=$#dp+1;

my %mqq;
my %fsq;
my %qdq;
my %qdq2;
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

$qdq{1}=$qd[sprintf("%.0f",$total*0.01)];
$qdq{5}=$qd[sprintf("%.0f",$total*0.05)];
$qdq{10}=$qd[sprintf("%.0f",$total*0.1)];
$qdq{15}=$qd[sprintf("%.0f",$total*0.15)];
$qdq{20}=$qd[sprintf("%.0f",$total*0.2)];
$qdq{30}=$qd[sprintf("%.0f",$total*0.3)];
$qdq{40}=$qd[sprintf("%.0f",$total*0.4)];
$qdq{50}=$qd[sprintf("%.0f",$total*0.5)];
$qdq{60}=$qd[sprintf("%.0f",$total*0.6)];
$qdq{70}=$qd[sprintf("%.0f",$total*0.7)];
$qdq{80}=$qd[sprintf("%.0f",$total*0.8)];
$qdq{90}=$qd[sprintf("%.0f",$total*0.9)];

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


print STDERR "\nMQ quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %mqq){
	print STDERR "$q\t$mqq{$q}\n";
}

print STDERR "\nFS quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %fsq){
	print STDERR "$q\t$fsq{$q}\n";
}

print STDERR "\nDP quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %dpq){
	print STDERR "$q\t$dpq{$q}\t$dpq2{$q}\n";
}

print STDERR "\nQD quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %qdq){
	print STDERR "$q\t$qdq{$q}\t$qdq2{$q}\n";
}

my @joint=();
my %jq={};

open TR, $true;
while(<TR>){
	if ($_=~/^#/) { 
		next;
	}
#warn "---------------\n";
	@dline=split("\t",$_);
	@iline=split(";",$dline[7]);
	my $mqQ=0;
	my $fsQ=0;
	my $dpQ=0;
	my $qdQ=0;
	my $mqt;
	my $fst;
	my $dpt;
	my $qdt;
	foreach my $i (@iline){
		if ($i=~/MQ=([-\.\d]+)/) { 
			if ($nomq){
				$mqQ=100;
				next;
			}
			$mqt=$1;
			foreach my $mqv (sort {$a <=> $b} keys %mqq) {
				if ($mqt<$mqq{$mqv}) { 
					$mqQ=$mqv;
					last;
				}
			}
			if (!$mqQ){$mqQ=100;}	
#warn "MQ=$mqt:$mqQ\n";		
		}
		elsif ($i=~/FS=([-\.\d]+)/) { 
			if ($nofs){
				$fsQ=100;
				next;
			}
			$fst=$1;
			foreach my $fsv (sort {$a <=> $b} keys %fsq) {
				if ($fst>$fsq{$fsv}) { 
					$fsQ=$fsv;
					last;
				}
			}			
			if (!$fsQ){$fsQ=100;}			
#warn "FS=$fst:$fsQ\n";		
		}
		elsif ($i=~/QD=([-\.\d]+)/) { 
			if ($noqd){
				$qdQ=100;
				next;
			}
			$qdt=$1;
			foreach my $qdv (sort {$a <=> $b} keys %qdq) {
				if ($qdt>$qdq{$qdv}) { 
					$qdQ=$qdv;
					last;
				}
			}			
			if (!$qdQ){$qdQ=100;}			
#warn "QD=$qdt:$qdQ\n";		
		}
		elsif ($i=~/DP=(\d+)/) { 
			if ($nodp){
				$dpQ=100;
				next;
			}
			$dpt=$1;
			foreach my $dpv (sort {$a <=> $b} keys %dpq) {
				if ($dpt<$dpq{$dpv} or $dpt>$dpq2{$dpv}) { 
					$dpQ=$dpv;
					last;
				}
			}
			if (!$dpQ){
				$dpQ=100;
			}
#			if ($dpQ>50) { $dpQ=100-$dpQ;}						
#warn "DP=$dpt:$dpQ\n";		
		}
	}
	push @joint, sprintf("%.5f",($dpQ/100)*($qdQ/100)*($mqQ/100)*($fsQ/100));
}
close TR;

@joint=sort {$a <=> $b} @joint;
$jq{1}=$joint[sprintf("%.0f",$total*0.01)];
$jq{5}=$joint[sprintf("%.0f",$total*0.05)];
$jq{10}=$joint[sprintf("%.0f",$total*0.1)];
$jq{15}=$joint[sprintf("%.0f",$total*0.15)];
$jq{20}=$joint[sprintf("%.0f",$total*0.2)];
$jq{30}=$joint[sprintf("%.0f",$total*0.3)];
$jq{40}=$joint[sprintf("%.0f",$total*0.4)];
$jq{50}=$joint[sprintf("%.0f",$total*0.5)];
$jq{60}=$joint[sprintf("%.0f",$total*0.6)];
$jq{70}=$joint[sprintf("%.0f",$total*0.7)];
$jq{80}=$joint[sprintf("%.0f",$total*0.8)];
$jq{90}=$joint[sprintf("%.0f",$total*0.9)];

print STDERR "\nJOINT quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %jq){
	next if ($q=~/HASH/);
	print STDERR "$q\t$jq{$q}\n";
}

my $ones=0;
my $fives=0;
my $tens=0;
my $fifteens=0;
my $twenties=0;
my $thirties=0;
my $total2=0;

open VCF, $vcf or die "cannot open vcf $vcf\n";

while(<VCF>){
	if ($_=~/^#/) { 
		print $_;
		next;
	}
	$total2++;
	@dline=split("\t",$_);
	@iline=split(";",$dline[7]);
	my $mqQ=0;
	my $fsQ=0;
	my $dpQ=0;
	my $jQ=0;
	foreach my $i (@iline){
		if ($i=~/MQ=(\d+)/) { 
			if ($nomq){
				$mqQ=100;
				next;
			}
			my $mqt=$1;
			foreach my $mqv (sort {$a <=> $b} keys %mqq) {
				if ($mqt<$mqq{$mqv}) { 
					$mqQ=$mqv;
					last;
				}
			}
			if (!$mqQ){$mqQ=100;}			
		}
		elsif ($i=~/FS=(\d+)/) { 
			if ($nofs){
				$fsQ=100;
				next;
			}
			my $fst=$1;
			foreach my $fsv (sort {$a <=> $b} keys %fsq) {
				if ($fst>$fsq{$fsv}) { 
					$fsQ=$fsv;
					last;
				}
			}			
			if (!$fsQ){$fsQ=100;}			
		}
		elsif ($i=~/QD=(\d+)/) { 
			if ($noqd){
				$qdQ=100;
				next;
			}
			my $qdt=$1;
			foreach my $qdv (sort {$a <=> $b} keys %qdq) {
				if ($qdt>$qdq{$qdv}) { 
					$qdQ=$qdv;
					last;
				}
			}			
			if (!$qdQ){$qdQ=100;}			
		}
		elsif ($i=~/DP=(\d+)/) { 
			if ($nodp){
				$dpQ=100;
				next;
			}
			my $dpt=$1;
			foreach my $dpv (sort {$a <=> $b} keys %dpq) {
				if ($dpt<$dpq{$dpv} or $dpt>$dpq2{$dpv}) { 
					$dpQ=$dpv;
					last;
				}
			}
			if (!$dpQ){
				$dpQ=100;
			}
#			if ($dpQ>50) { $dpQ=100-$dpQ;}						
		}
	}
	my $jt=sprintf("%.5f",($dpQ/100)*($qdQ/100)*($mqQ/100)*($fsQ/100));
#warn "\t$jt\n";
	foreach my $jv (sort {$a <=> $b} keys %jq) {
		if ($jt<$jq{$jv}) { 
			$jQ=$jv;
			last;
		}
	}			
	if (!$jQ){$jQ=100;}
	if ($jQ==1) { $ones++;}
	elsif ($jQ==5) { $fives++;}
	elsif ($jQ==10) { $tens++;}			
	elsif ($jQ==15) { $fifteens++;}			
	elsif ($jQ==20) { $twenties++;}			
	elsif ($jQ==30) { $thirties++;}	
	$jQ=$jQ-0.1;		
	$dline[5]=$jQ;
	push @iline,"SBQ=".$fsQ;
	push @iline,"ABQ=".$mqQ;
	push @iline,"DPQ=".$dpQ;
	push @iline,"TPQ=".$qdQ;
	$dline[7]=join(";",@iline);	
	print join("\t",@dline);
}

my $oness=sprintf("%.2f",100*$ones/$total2);
my $fivess=sprintf("%.2f",100*($ones+$fives)/$total2);
my $tenss=sprintf("%.2f",100*($ones+$fives+$tens)/$total2);
my $fifteenss=sprintf("%.2f",100*($ones+$fives+$tens+$fifteens)/$total2);
my $twentiess=sprintf("%.2f",100*($ones+$fives+$tens+$fifteens+$twenties)/$total2);
my $thirtiess=sprintf("%.2f",100*($ones+$fives+$tens+$fifteens+$twenties+$thirties)/$total2);

print STDERR "
------------------------
$oness%\tat qual <1 (",sprintf("%.2f",$oness-1),"% gain)
$fivess%\tat qual <5 (",sprintf("%.2f",$fivess-5),"% gain)
$tenss%\tat qual <10 (",sprintf("%.2f",$tenss-10),"% gain)
$fifteenss%\tat qual <15 (",sprintf("%.2f",$fifteenss-15),"% gain)
$twentiess%\tat qual <20 (",sprintf("%.2f",$twentiess-20),"% gain)
$thirtiess%\tat qual <30 (",sprintf("%.2f",$thirtiess-30),"% gain)
------------------------

";


