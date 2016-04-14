#!/usr/bin/perl -w
use lib "/home/dhendrix/PerlModules/";
use miRTRAP;
use strict;
$|=1;

my $parameters = {};
my $usage = "Usage:\n$0 <config file>\n";
my $configFile = $ARGV[0] or die $usage;
miRTRAP::readConfigFile($configFile,$parameters);

my $scriptDir = $parameters->{scriptDir} ? $parameters->{scriptDir} : "";
my $printReadRegions = $scriptDir . "/" . "printReadRegions.pl";
my $readRegionsToRNAFold = $scriptDir . "/" . "readRegionsToRNAFold.pl";
my $RNAfold = "RNAfold";
my $trnascan = "tRNAscan-SE";
my $processReadRegions = $scriptDir . "/" ."processReadRegions.pl";
my $evaluateReadRegionsWithGenes = $scriptDir . "/" . "evaluateReadRegionsWithGenes.pl";
my $filePrefix = $parameters->{filePrefix} ? $parameters->{filePrefix} : "readRegions";

runPrintReadRegions($parameters);
runReadRegionsToRNAFold($parameters);
runRNAFold($parameters);
runTRNASCAN($parameters);
runProcessReadRegions($parameters);
evaluateReadRegionsWithGenes($parameters);

sub runPrintReadRegions {
    my($parameters) = @_;
    my $readListFile = $parameters->{readListFile} or die "FAIL: readListFile not loaded in config file.\n";
    my $repeatRegionsFile = $parameters->{repeatRegionsFile} or die "FAIL: repeatRegionsFile not loaded in config file.\n";
    my $configFile = $parameters->{configFile} or die "FAIL: configFile not loaded.\n";
    open(RLF,$readListFile) or die "FAIL: could not open $readListFile\n";
    close(RLF);
    # repeat regions are not required.
    #open(RRF,$repeatRegionsFile) or die "FAIL: could not open $repeatRegionsFile\n";
    #close(RRF);
    system("$printReadRegions $configFile");
}

sub runReadRegionsToRNAFold {
    my($parameters) = @_;
    my $readRegionsFile = $parameters->{readRegionsFile} or die "FAIL: readRegionsFile not found.\n";
    my $genomeFile = $parameters->{genomeFile} or die "FAIL: genomeFile not loaded in config file.\n";    
    my $readRegionsFastaFile = $parameters->{readRegionsFastaFile} or die "FAIL: readRegionsFastaFile not found.\n";
    my $configFile = $parameters->{configFile} or die "FAIL: configFile not loaded.\n";
    open(RRF,$readRegionsFile) or die "FAIL: could not open $readRegionsFile\n";
    close(RRF);
    open(GF,$genomeFile) or die "FAIL: could not open $genomeFile\n";
    close(GF);
    system("$readRegionsToRNAFold $configFile");
}

sub runRNAFold {
    my($parameters) = @_;
    my $readRegionsFastaFile = $parameters->{readRegionsFastaFile} or die "FAIL: readRegionsFastaFile not found.\n";
    my $rnaFoldOutput = $filePrefix.".mfe";
    open(RRF,$readRegionsFastaFile) or die "FAIL: could not open $readRegionsFastaFile\n";
    close(RRF);
    system("cat $readRegionsFastaFile | $RNAfold -noPS > $rnaFoldOutput");
    return $rnaFoldOutput;
}

sub runTRNASCAN {
    my($parameters) = @_;
    my $readRegionsFastaFile = $parameters->{readRegionsFastaFile} or die "FAIL: readRegionsFastaFile not found.\n";
    my $trnascanOutputFile = $filePrefix.".trna";
    open(RRF,$readRegionsFastaFile) or die "FAIL: could not open $readRegionsFastaFile\n";
    close(RRF);
    system("$trnascan -o $trnascanOutputFile -q $readRegionsFastaFile");
    $parameters->{trnascanOutputFile} = $trnascanOutputFile;
}

sub runProcessReadRegions {
    my($parameters) = @_;
    my $readListFile = $parameters->{readListFile} or die "FAIL: readListFile not loaded in config file.\n";
    my $rnaFoldOutput = $parameters->{rnaFoldOutput} or die "FAIL: cound not load rnaFoldOutput.\n";
    my $trnascanOutputFile = $parameters->{trnascanOutputFile} or die "FAIL: trnascanOutputFile not found.\n";    
    my $configFile = $parameters->{configFile} or die "FAIL: configFile not loaded.\n";
    open(RNAF,$rnaFoldOutput) or die "FAIL: could not open $rnaFoldOutput\n";
    close(RNAF);
    open(TRNAS,$trnascanOutputFile) or die "FAIL: could not open $trnascanOutputFile\n";
    close(TRNAS);
    system("$processReadRegions $configFile");
}

sub evaluateReadRegionsWithGenes {
    my($parameters) = @_;
    my $hairpinListFile = $parameters->{hairpinsFile} or die "hairpinListFile: parameter not loaded.\n";
    my $hairpinProductsFile = $parameters->{productFile} or die "hairpinProductsfile: parameter not loaded.\n";
    my $readRegionsFile = $parameters->{allReadRegionsFile} or die "allReadRegionsFile: parameter not loaded.\n";
    my $configFile = $parameters->{configFile} or die "FAIL: configFile not loaded.\n";
    my $geneModels = $parameters->{geneModels} or die "geneModels: parameter not loaded.\n";
    system("$evaluateReadRegionsWithGenes $configFile $geneModels");
}
