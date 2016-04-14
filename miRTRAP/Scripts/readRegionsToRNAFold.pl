#!/usr/bin/perl -w
use miRTRAP;
use strict;
$|=1;

my $usage = "Usage:\n\n$0 <config file>\n";

my $parameters = {
    "totalLength" => 150,
    "filePrefix" => "readRegions"
};

my $configFile = $ARGV[0] or die $usage;
miRTRAP::readConfigFile($configFile,$parameters);
my $readRegionsFile = $parameters->{readRegionsFile} or die "FAIL: readRegionsFile not loaded from configFile.\n";
my $genomeFile = $parameters->{genomeFile} or die "FAIL: genomeFile not loaded from configFile.";
my $chromRegExp = $parameters->{chromRegExp};
my $genome = miRTRAP::loadGenome($genomeFile,$chromRegExp);
miRTRAP::printReadRegionRNAfold($readRegionsFile,$genome,$parameters);
