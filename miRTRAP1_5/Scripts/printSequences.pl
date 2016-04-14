#!/usr/bin/perl -w
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 <fasta or fastq>\n";
my $sequencesFile = $ARGV[0] or die $USAGE;

if ($sequencesFile =~ /fastq$/) {
    fastqToCondensedFasta($sequencesFile);
} elsif ($sequencesFile =~ /fasta$/) {
    fastaToCondensedFasta($sequencesFile);
} else {
    die "$sequencesFile has unknown file type\n";
}

sub fastqToCondensedFasta {
    my($fastqFile) = @_;
    open(FQ,$fastqFile);
    while (<FQ>) {
	my $def = $_;
	my $seq = <FQ>;
	my $def2 = <FQ>;
	my $qual = <FQ>;
	print $seq unless (length($seq) < 15);
    }
}

sub fastaToCondensedFasta {
    my($fastaFile) = @_;
    open(FA,$fastaFile);
    while (<FA>) {
	my $def = $_;
	my $seq = <FA>;
	print $seq unless (length($seq) < 15);
    }
}
