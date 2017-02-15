#!/bin/perl
############################################################
#Author: Gokul
#perl script that runs blat for all aligned reads for each editing sites
#inputs: List of variants, INDEXED!! BAM alignment file, Output file name

use warnings;
use strict;
use Data::Dumper;

if (@ARGV != 3) {
	die "need to provide 3 input:Variant list, INDEXED BAM alignment file and output file name\n";
}
my ($inputfile, $bamfile, $outputfile) = ($ARGV[0], $ARGV[1], $ARGV[2]);

my $minbasequal = 25;
my $minmismatch = 2;
my $scorelimit = 0.95;


open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);

my $fafile = join '', $inputfile, '.fatmp';
my $pslfile = join '', $inputfile, '.psltmp';
open (FAFILE, ">", $fafile);

while (<$INPUT>) {
	chomp;
	my $inputline = $_;
	my @fields = split;
	my $TEMPNAME = join '', $outputfile,'_tmp';
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $bamposition = join '', $chr,':',$position,'-',$position;
	system("~/liucj/tools/samtools-1.3.1/install/bin/samtools view $bamfile $bamposition > $TEMPNAME");
	my $editnuc = $fields[4];
	my $newmismatch = 0;
	my $mismatchreadcount = 0;

	open(my $TMPFILE, "<", $TEMPNAME);
	while (<$TMPFILE>) {
		chomp;
		my @bamfields = split;
		my ($alignment, $readstart, $cigar, $sequence, $qualities) = ($bamfields[1], $bamfields[3], $bamfields[5], $bamfields[9], $bamfields[10]);
		my @sequencebases = split(//,$sequence);
		my @qualscores = split(//,$qualities);
		my ($currentpos, $readpos) = ($readstart,1);
		my $base_readpos;
		my @cigarnums = split(/[MIDNSHP]/, $cigar);
		my @cigarletters = split(/[0-9]+/,$cigar);
		shift @cigarletters;
		for (my $i = 0; $i < @cigarnums; $i++) {
			if ($cigarletters[$i] =~ m/[ISH]/) {
				$readpos = $readpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/[DN]/) {
				$currentpos = $currentpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/M/) {
				for (my $j = 0; $j < $cigarnums[$i]; $j++) {
					$base_readpos = 1 if ($currentpos == $position && $sequencebases[$readpos-1] eq $editnuc && ord($qualscores[$readpos-1]) >= $minbasequal+33);
					$currentpos++;
					$readpos++;	
				}	
			}
		}
		if ($base_readpos) {
			print FAFILE ">$chr\-$position\-$mismatchreadcount\n$sequence\n";
			$mismatchreadcount++;
		}
	}
	close $TMPFILE;
}
 close $INPUT;


# system("/home/cliu18/liucj/tools/kent/src/blat/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -noHead /home/cliu18/liucj/pipelines/exome_pipeline/data/hg38/genomeBuild/hg38.fasta $fafile $pslfile");
# system("rm $fafile");




open(PSL, "<", $pslfile );
my %pslhash;
while(<PSL>) {
	chomp;
	my @pslfields = split;
	my $name = $pslfields[9];
	my $blatscore = join '@',$pslfields[0],$pslfields[13],$pslfields[17],$pslfields[18],$pslfields[20];
	if ($pslhash{$name}) {
		$pslhash{$name} = join '-', $pslhash{$name}, $blatscore;
	} elsif ( !$pslhash{$name}) {
		$pslhash{$name} = $blatscore; 
	}
}
close PSL;	
print Dumper(\%pslhash);

my %sitehash;
my %discardhash;
foreach my $pslkey (keys %pslhash) {
	my @splitkey = split(/\-/, $pslkey);
	my $site = join '_', $splitkey[0],$splitkey[1];
	my @psllines = split(/\-/, $pslhash{$pslkey});
	my $largestscore = 0;
	my $largestscoreline = $psllines[0];
	my @scorearray;
	foreach my $scoreline (@psllines) {
		my @scoresarray = split(/\@/,$scoreline);
		my $linescore = $scoresarray[0];
		push(@scorearray,$linescore);
		if ($linescore > $largestscore) {
			$largestscoreline = $scoreline;
			$largestscore = $linescore;
		}
	}
	@scorearray = sort {$b <=> $a} @scorearray;
	$scorearray[1] = 0 unless ($scorearray[1]);
    # print join "\t", @scorearray, "\n";
    # print $largestscoreline, "\n";
	my @splitlargestline = split(/\@/,$largestscoreline);
	my $overlapfound = 0;
	if ($splitlargestline[1] eq $splitkey[0] && $scorearray[1] < ($scorearray[0]*$scorelimit)) {
		my ($numblocks, $blocksizes, $blockstarts) = ($splitlargestline[2],$splitlargestline[3],$splitlargestline[4]);
		my @blocksizes = split(/\,/,$blocksizes);
		my @blockstarts = split(/\,/,$blockstarts);
		for (my $i = 0; $i < $numblocks; $i++) {
			my $startpos = $blockstarts[$i]+1;
			my $endpos = $blockstarts[$i] + $blocksizes[$i];
			$overlapfound = 1 if ($splitkey[1] >= $startpos && $splitkey[1] <= $endpos);
		}
		if ($sitehash{$site} && $overlapfound) {
			$sitehash{$site}++;
		} elsif (!$sitehash{$site} && $overlapfound) {
			$sitehash{$site} = 1;
		}
	}
	unless ($overlapfound) {
		if ($discardhash{$site}) {
			$discardhash{$site}++;
		} elsif (!$discardhash{$site}) {
			$discardhash{$site} = 1;
		}
	}
}
print Dumper(\%sitehash);
print Dumper(\%discardhash);

open (SITES2, "<", $inputfile ) or die "error opening inputfile: $!\n";
while(<SITES2>) {
	chomp;
	my @fields = split;
	my $inputline = $_;
	my ($cov,$oldalter) = split(/\,/,$fields[2]);
	my $name = join '', $fields[0],'_',$fields[1];
	if ($sitehash{$name}) {
		my $newalter = $sitehash{$name};
		my $discardnum;
		if ($discardhash{$name}) {
			$discardnum = $discardhash{$name};
		} else {
			$discardnum = 0;
		}
		my $newcov = $cov - ($oldalter - $newalter);
		my $neweditfreq = sprintf("%.3f", $newalter/$newcov);
		print $OUTPUT "$fields[0]\t$fields[1]\t$newcov,$newalter\t$fields[3]\t$fields[4]\t$neweditfreq\t$fields[6]\t$fields[7]\t$fields[8]\n" if ($newalter >= $minmismatch && $newalter > $discardnum);
		#print $OUTPUT "$inputline\n" if ($discardnum < $minmismatch && $newalter >= $minmismatch);
	}
}
close SITES2;
close $OUTPUT;
# system("rm $pslfile");
