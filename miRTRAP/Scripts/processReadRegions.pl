#!/usr/bin/perl -w
use miRTRAP;
use strict;
$|=1;

my $usage = "Usage:\n\n$0 <configFile>\n";

my $parameters = {
    "minLength" => 20,
    "minDist" => 10,
    "minMajor" => -44,
    "maxMajor" => 22,
    "maxReverse" => 0.05,
    "maxHitCount" => 50,
    "minLocusCount" => 5,
    "maxFivePrimeHet" => 0.5,
    "minShift" => 7,
    "minOverlap" => 2,
    "bpDensityLimit" => 0.6,
    "inHairpinBuffer" => 3,
    "outHairpinBuffer" => 3,
    "hairpinRange" => 70,
    "RNAfold" => "RNAfold"
    };

my $configFile = $ARGV[0] or die $usage;
miRTRAP::readConfigFile($configFile,$parameters);
my $rnaFoldOutput = $parameters->{rnaFoldOutput} or die "rnaFoldOutput file not found.\n";
my $readGffListFile = $parameters->{readListFile} or die "readListFile file not found.\n";
my $tRNAScanFile = $parameters->{trnascanOutputFile} ? $parameters->{trnascanOutputFile} : "";
my($dataSet,$sampleTotal,$sampleList,$hitCount) = miRTRAP::loadReadGffList($readGffListFile,$parameters);
my $foldData = miRTRAP::readRNAFoldOutput($rnaFoldOutput);
my $tRNAs = miRTRAP::readtRNAScanFile($tRNAScanFile);
miRTRAP::processReadRegions($foldData,$dataSet,$sampleList,$hitCount,$tRNAs,$parameters);
