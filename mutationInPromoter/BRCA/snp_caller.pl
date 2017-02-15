
#script used to call SNPs from RNA-SEQ pileup data

use warnings;
use strict;
require "/home/cliu18/liucj/github/RstudioWithGit/mutationInPromoter/BRCA/parse_pileup.pl";

if (@ARGV != 1) {
	die "need to provide 1 input:pileupfile\n";
}
my ($mincov, $maxcov, $minvarfreq, $minbasequal, $minvarcount) = (1, 15000000, 0.1, 20, 2); #set variant calling parameters
my $pileupfile = $ARGV[0];
open (PILEUP, "gunzip < $pileupfile |") or die "error opening pileup input: $!\n"; #unzip pileup file line by line and stream into program
while(<PILEUP>) {
	chomp;
	my ($chr, $position, $refnuc, $coverage, $pile, $qual) = split;
	my ($refnuccount, $acount, $tcount, $ccount, $gcount) = &parse_pileup($_, $minbasequal);# parse each line of pileup
	$coverage = $refnuccount + $acount + $tcount + $ccount + $gcount + 0.0; 
	if ($coverage >= $mincov && $coverage <= $maxcov) { #check coverage is within bounds
		my $snp = 0;
		my $mutnuc;
		my $varfreq;
		my $mutnuccount;
		if ($acount >= $minvarcount && $acount/$coverage >= $minvarfreq) { #check for variant frequency >= 10% and variant count >=2 for each nucleotide
			$snp = 1;
			$mutnuc = 'A';
			$varfreq = sprintf("%.3f", $acount/$coverage);
			$mutnuccount = $acount;
		}
		elsif ($tcount >= $minvarcount && $tcount/$coverage >= $minvarfreq) {
			$snp = 1;
			$mutnuc = 'T';
			$varfreq = sprintf("%.3f", $tcount/$coverage);
			$mutnuccount = $tcount;
		}
		elsif ($ccount >= $minvarcount && $ccount/$coverage >= $minvarfreq) {
			$snp = 1;
			$mutnuc = 'C';
			$varfreq = sprintf("%.3f", $ccount/$coverage);
			$mutnuccount = $ccount;
		}
		elsif ($gcount >= $minvarcount && $gcount/$coverage >= $minvarfreq) {
			$snp = 1;
			$mutnuc = 'G';
			$varfreq = sprintf("%.3f", $gcount/$coverage);
			$mutnuccount = $gcount;
		}
		print "$chr\t$position\t$coverage,$mutnuccount\t$refnuc\t$mutnuc\t$varfreq\n" if ($snp); # if SNP exists print out info
	}

}
close PILEUP;
