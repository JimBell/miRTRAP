#!/usr/bin/perl -w
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 <fasta or fastq file>\n";

my $sequencesFile = $ARGV[0] or die $USAGE;
my($fileBase) = $sequencesFile =~ /([^\/]*)\.fast.$/;
my $outputFile = $fileBase . "_condensed.fasta";

convertFastqToCondensedFasta($sequencesFile,$outputFile);


sub convertFastqToCondensedFasta {
    my($sequencesFile,$outputFile) = @_;
    open(SORTEDSEQS, "~/Scripts/miRTRAP1_5/printSequences.pl $sequencesFile | sort | uniq -c |");
    open(FASTA,">$outputFile");
    my $readNum = 1;
    while (<SORTEDSEQS>) {
	chomp;
	my($count,$seq) = split;
	print FASTA ">dr".$readNum."_x".$count."\n";
	print FASTA "$seq\n";
	$readNum++;
    }
    close(SORTEDSEQS);
    close(FASTA);
}


