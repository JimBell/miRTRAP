#!/usr/bin/perl -w
use miRTRAP;
use strict;
$|=1;

my $usage = "Usage:\n\n$0 <configFile>\n";

my $parameters = {
    "maxLength" => 160,
    "maxCount" => 5,       
    "maxHitCount" => 50,
    "filePrefix" => "readRegions",
    "readListFile" => "",
    "repeatRegionsFile" => ""	
    };

my $configFile = $ARGV[0] or die $usage;
miRTRAP::readConfigFile($configFile,$parameters);
my $readGffListFile = $parameters->{readListFile} or die "FAIL: readListFile not loaded in configFile.\n";
my $repeatRegionsFile = $parameters->{repeatRegionsFile};
my($dataSet,$sampleTotal,$sampleList,$hitCount) = miRTRAP::loadReadGffList($readGffListFile,$parameters);
my $repeatRegions = {};
if($repeatRegionsFile) {
    $repeatRegions = miRTRAP::readRepeatRegions($repeatRegionsFile);
}
miRTRAP::processOverlaps($dataSet,$repeatRegions,$hitCount,$parameters);
