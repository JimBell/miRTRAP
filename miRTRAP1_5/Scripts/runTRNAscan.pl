#!/usr/bin/perl
use miRTRAP1_5;
use Getopt::Long;
use strict;

$| = 1;

unless (@ARGV) {
    miRTRAP1_5::printUsage();
    die;
}

my $parameters = miRTRAP1_5::loadDefaultParameters();

Getopt::Long::Configure("no_ignore_case");

GetOptions ("LengthMax=i" => \$parameters->{LengthMax},
	    "lengthMin=i" => \$parameters->{lengthMin},
	    "totalLength=i" => \$parameters->{totalLength},
	    "distanceMin=i" => \$parameters->{distanceMin},
	    "CountMax=i" => \$parameters->{CountMax},
	    "HitCountmax=i" => \$parameters->{HitCountMax},
	    "MaxMajor=i" => \$parameters->{MaxMajor},
	    "minMajor=i" => \$parameters->{minMajor},
	    "reverseMax=i" => \$parameters->{reverseMax},
	    "countMinLocus=i" => \$parameters->{countMinLocus},
	    "fivePrimeHetMax=f" => \$parameters->{fivePrimeHetMax},
	    "shiftMin=i" => \$parameters->{shiftMin},
	    "OverlapMin=i" => \$parameters->{OverlapMin},
	    "BPdensityLimit=f" => \$parameters->{BPdensityLimit},
	    "HairpinBufferIn=i" => \$parameters->{HairpinBufferIn},
	    "hairpinBufferOut=i" => \$parameters->{hairpinBufferOut},
	    "RangeOfHairpin=i" => \$parameters->{RangeOfHairpin},
	    "outputPrefix=s" => \$parameters->{outputPrefix},
	    "bamListFile=s" => \$parameters->{bamListFile},
	    "FileRepeatRegions=s" => \$parameters->{FileRepeatRegions},
	    "genomeDir=s" => \$parameters->{genomeDir},
	    "SizesFile=s" => \$parameters->{SizesFile}
    );

miRTRAP1_5::createOutputFileParameters($parameters);
my $tRNAScanFastaFile = $parameters->{tRNAScanFasta};
my $tRNAScanOutputFile = $parameters->{tRNAScanOutputFile};
my $hairpinsFile = $parameters->{hairpinsFile};
miRTRAP1_5::printTRNAFastaFile($hairpinsFile, $tRNAScanFasta);
System("tRNAScan-SE -o $tRNAScanOutputFile -q $hairpinsFile");

