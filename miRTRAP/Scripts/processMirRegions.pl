#!/usr/bin/perl -w
use miRTRAP;
use strict;
$|=1;

my $parameters = {
    "minLength" => 10,
    "minDist" => 10,
    "minMajor" => -44,
    "maxMajor" => 22,
    "maxReverse" => 0.05,
    "maxHitCount" => 50,
    "maxFivePrimeHet" => 0.5,
    "minShift" => 7,
    "minOverlap" => 2,
    "bpDensityLimit" => 0.5,
    "inHairpinBuffer" => 3,
    "outHairpinBuffer" => 10,
    "hairpinRange" => 70,
    "RNAfold" => "RNAfold"
    };

my $usage = "Usage:\n\n$0 <mir regions file> <deep seq read GFF list>\n";
my $mirRegionsFile = $ARGV[0] or die $usage;
my $readGffListFile = $ARGV[1] or die $usage;
my($dataSet,$sampleTotal,$sampleList,$hitCount) = miRTRAP::loadReadGffList($readGffListFile,$parameters);
my $foldData = miRTRAP::readMirRegionsFile($mirRegionsFile);
miRTRAP::processMirRegions($foldData,$dataSet,$sampleList,$hitCount,$parameters);
