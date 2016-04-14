#!/usr/bin/perl -w
use miRTRAP1_6;
use Bio::DB::Sam;
use Getopt::Long;
use strict;
$|=1;

my $USAGE = "USAGE:\n$0 -S <SizesFile> -b <BamListFile> -o <outputPrefix> [-R <Repeat regions file>]
\t-S\tSizesFile\tchrom.sizes file containing the size of each chromosome.
\t-b\tbamListFile\t\tfile containing a list of sample names and bamfile locations
\t\t\t\tbamListFileExample:
\t\t\t\t\t<sample name 1>\t<bam file 1>
\t\t\t\t\t<sample name 2>\t<bam file 2>
\t\t\t\t\t<sample name 3>\t<bam file 3>
\t\t\t\t\tetc.
\t-R\tRepeatRegionsFile\tFile containing a list of repeat regions. (this file is optional)
\t-o\toutputPrefix\tprefix for output files;
\t-M\tMaxLength\tThe maximum length for each read region reported in the read regions file\n";

die $USAGE unless (@ARGV);

my $parameters = miRTRAP1_6::loadDefaultParameters();

Getopt::Long::Configure("no_ignore_case");

GetOptions ("MaxLength=i" => \$parameters->{MaxLength},
	    "outputPrefix=s" => \$parameters->{outputPrefix},
	    "bamListFile=s" => \$parameters->{bamListFile},
	    "RepeatRegions=s" => \$parameters->{RepeatRegionsFile},
	    "SizesFile=s" => \$parameters->{SizesFile},
	    "LoadFromConfigFile=s" => \$parameters->{LoadFromConfigFile}
    );

unless ($parameters->{LoadFromConfigFile} eq "") {
    my $configFile = $parameters->{LoadFromConfigFile};
    miRTRAP1_6::readConfigFile($configFile,$parameters);
}
miRTRAP1_6::createOutputFileParameters($parameters);

my $bamListFile = $parameters->{bamListFile} or die "FAIL: Bam List file not loaded (not found in parameters).\n";
my $chromSizesFile = $parameters->{SizesFile} or die "FAIL: chrom.sizes file not loaded (not found in parameters).\n";
my $repeatRegionsFile = $parameters->{RepeatRegionsFile};
my $chromLengths = miRTRAP1_6::readChromLengths($chromSizesFile) or die "Failed to load chrom.sizes file.\n";
my $bamList = miRTRAP1_6::loadBamList($bamListFile) or die "Failed to load bamList file.\n";
my $repeatRegions = {};
if($repeatRegionsFile) {
    $repeatRegions = miRTRAP1_6::loadRepeatRegions($repeatRegionsFile, $chromLengths) or die "Failed to load repeat regions file.\n";
}
miRTRAP1_6::printReadRegions($bamList, $chromLengths, $repeatRegions, $parameters);
