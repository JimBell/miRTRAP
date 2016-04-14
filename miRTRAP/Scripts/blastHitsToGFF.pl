#!/usr/bin/perl -w
use POSIX qw(floor);
use strict;
$|=1;

my $usage = "$0 <blast hits file list> <source>\n";

my $fileList = $ARGV[0] or die $usage;
my $source = $ARGV[1] or die $usage;
my $type = "match";
my $eValueLimit=0.01;

open(FL,$fileList) or die "could not open file: $fileList\n";
open(AF,">".$source."_reads.gff");
while(<FL>) {
    chomp;
    open(BHF,$_) or die "could not open file: $_\n";
    while(<BHF>) {
	chomp;
	my($qCoord,$hits,$sCoord,$eValue,$score,$length,$pid)=split(/ /,$_);
	if($eValue <= $eValueLimit) {
	    my($sRef,$sStart,$sStop,$sStrand)=parseLocation($sCoord);
	    my($qRef,$qStart,$qStop,$qStrand)=parseLocation($qCoord);
	    #$sRef =~ s/chr//g;
	    $sRef =~ s/scaffold/Scaffold/g;
	    print AF "$sRef\t$source\t$type\t$sStart\t$sStop\t$eValue\t$sStrand\t.\tMatch $qRef\n";
	}
    }
    close(BHF);
}
close(AF);

sub parseLocation {
    my($location)=@_;
    my($chrom,$start,$end,$strand);
    if($location =~ /(.*)\:(\d+)\-(\d+)\:(.*)/) {	
	$chrom=$1;
	$start=$2;
	$end=$3;
	$strand=$4;
    } elsif($location =~ /(.*)\:(\d+)\-(\d+)/) {
	$chrom=$1;
	$start=$2;
	$end=$3;
	$strand="+";
    }
    return ($chrom,$start,$end,$strand);
}
