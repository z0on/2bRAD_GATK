#!/usr/bin/perl

$usage= "

ngs_concat.pl : 
concatenates files by matching pattern in their names

arg1: common pattern for files
arg2: perl-like pattern in the filename to recognize, 
	  use brackets to specify the unique part

example: ngs_concat.pl 'Sample' 'Sample_(..)'

";

my $ff = shift or die $usage;
my $patt=shift or die $usage;
#print "pattern $patt\n";
opendir THIS, ".";
my @files=grep /$ff/,readdir THIS;
print "files:\n",@files, "\n";
my @ids=();
foreach $file (@files){
	if ($file=~/$patt/) {
		$ii=$1;
#		unless (grep {$_ eq $ii} @ids){ 
#			my $name=$file;
#			my $ccat=$patt;
#			$ccat=~s/\(.+\)/$ii/;
#			$ccat.="*";
			$name=$ii.".fq";
print "$file > $name\n";

			`cat $file >> $name`;
#			push @ids, $ii;
#		}
	}
	else { print "$patt not found in $file\n";}
}

