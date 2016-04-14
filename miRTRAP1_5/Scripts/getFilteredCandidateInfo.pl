#!/usr/bin/perl -w
use Memory::Usage;
use miRTRAP1_5;
use Bio::DB::Sam;
use Getopt::Long;
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 -S <SizesFile> -b <BamListFile> -g <genome directory> -o <outputPrefix>
\t-S\tSizesFile\tchrom.sizes file containing the size of each chromosome.
\t-b\tbamListFile\t\tfile containing a list of sample names and bamfile locations
\t\t\t\tbamListFileExample:
\t\t\t\t\t<sample name 1>\t<bam file 1>
\t\t\t\t\t<sample name 2>\t<bam file 2>
\t\t\t\t\t<sample name 3>\t<bam file 3>
\t\t\t\t\tetc.
\t-g\tgenomeDir\tDirectory containing chromosomes in sperate fasta files
\t-o\toutputPrefix\tprefix for output files;
\t-l\tlengthMin\tminimum length for the arms of the hairpin
\t-t\ttotalLength\tmaximum length of the entire hairpin
\t-d\tdistanceMin\tdistance within which a read can still be associated with a product
\t-f\tfivePrimeHetMax\tMaximum 5` hererogeneity
\t-c\tcountMinLocus\tmaximum number of producs reads neccessary to be a real product (use total of values returned from addNHTags.pl)
\t-r\treverseMax\tmaximum allowed fraction of total antisense product reads to total (sense + antisense) product reads
\t-s\tshiftMin\tminimum shift of products
\t-h\thairpinShortLength\tlength of short arms allowed in the middle of hammer head loops
\t-O\tOverlapMin\tmin amount of product overlap needed for a same shifted and both shifted value to be recorded
\t-I\tInHairpinBuffer\tamount 5p products are allowed to cross into the loop without being considered a loop product
\t-O\tOutHairpinBuffer\tamount 3p products are allowed to cross into the loop without being considered a loop product
\t-R\tRangeOfHairpint\tRange of hairpin arms.  Outside this range are out products.
\t-M\tMaxLength\tThe maximum length for each read region reported in the read regions file\n";


die $USAGE unless (@ARGV);

my $parameters = miRTRAP1_5::loadDefaultParameters();

Getopt::Long::Configure("no_ignore_case");

GetOptions ("lengthMin=i" => \$parameters->{lengthMin},
	    "totalLength=i" => \$parameters->{totalLength},
	    "hairpinShortLength=i" => \$parameters->{hairpinShortLength},
	    "distanceMin=i" => \$parameters->{distanceMin},
	    "reverseMax=i" => \$parameters->{reverseMax},
	    "countMinLocus=i" => \$parameters->{countMinLocus},
	    "fivePrimeHetMax=f" => \$parameters->{fivePrimeHetMax},
	    "shiftMin=i" => \$parameters->{shiftMin},
	    "OverlapMin=i" => \$parameters->{OverlapMin},
	    "InHairpinBuffer=i" => \$parameters->{InHairpinBuffer},
	    "OutHairpinBuffer=i" => \$parameters->{OutHairpinBuffer},
	    "RangeOfHairpin=i" => \$parameters->{RangeOfHairpin},
	    "outputPrefix=s" => \$parameters->{outputPrefix},
	    "bamListFile=s" => \$parameters->{bamListFile},
	    "FileRepeatRegions=s" => \$parameters->{FileRepeatRegions},
	    "genomeDir=s" => \$parameters->{genomeDir},
	    "SizesFile=s" => \$parameters->{SizesFile}
    );

miRTRAP1_5::createOutputFileParameters($parameters);


my $bamListFile = $parameters->{bamListFile};
my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
my $chromSizesFile = $parameters->{SizesFile};
my $genomeDir = $parameters->{genomeDir};
my $chromLengths = miRTRAP1_5::readChromLengths($chromSizesFile);
my $bamList = miRTRAP1_5::loadBamList($bamListFile);
getFilteredCandidatesInfo($bamList,$filteredCandidatesFile,$genomeDir,$chromLengths,$parameters);

sub getFilteredCandidatesInfo {
    my($bamList,$filteredCandidatesFile,$genomeDir,$chromLengths,$parameters) = @_;
    my @sampleList;
    my $outputPrefix = $parameters->{outputPrefix};
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $hairpinsFile = $outputPrefix . "_filteredHairpins.txt";
    my $distinctReadsFile = $outputPrefix . "_filteredDistinctreads.txt";
    my $productFile = $outputPrefix . "_filteredProducts.txt";
    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
    open(HPL,">".$productFile) or die "failed to open $productFile for writing";
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tafh\tvalid\tsameShifted\tbothShift\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($sample) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";

    open(FRC, $filteredCandidatesFile) or die "failed to open $readRegionsFile\n";
    my $prevChrom;
    my $genome;
    while (<FRC>) {
        chomp;
	unless(/\#/) {	  
	    my($mirName,$location,$rejectionReason) = split(/\t/,$_);
	    my($chrom,$start,$stop,$strand) = miRTRAP1_5::parseLocation($location);
	    #the chroms in RRF should be sorted at this point so the next line of code will only be executed once per chrom
	    if ($prevChrom) {
		$genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	    } else {
		$genome = loadGenome("$genomeDir/$chrom.fa");	
	    }
	    $prevChrom = $chrom;
	    my $sequence = getSequence($location,$genome);
	    my($fold,$mfe) = RNA::fold($sequence);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList);

	    my @centers = getMergedHairpinCenters($fold,$parameters);
	    my $COUNT= 0;
	    foreach my $center (@centers) {
		my $asciiVal = ord('a') + $COUNT;
		my $label = chr($asciiVal);
		my $newId = $id . $label;
		$COUNT++;
		my $basePairs = getBasePairs($center,$fold);
		my $hairpinLength = getMaxHairpinLength($fold,$center);
		my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		if ($distinctReads->{$strand}) {
		    my $productInfo=extractProducts($center,
						    $fold,$basePairs,$location,
						    $distinctReads->{$strand},$strand,$parameters);
		    my $revStrand = revStrand($strand);
		    my $revProductInfo=extractProducts($center,
						       $fold,$basePairs,$location,
						       $distinctReads->{$revStrand},$revStrand,$parameters);
		    my $adjTotalProductReads = getAdjTotalProductReads($productInfo);
		    my $adjTotalRevProductReads = getAdjTotalProductReads($revProductInfo);
		    my($tpd,$totalRP)=getReverseProductDisplacement($productInfo,
								    $revProductInfo,
								    $adjTotalProductReads,
								    $parameters);
		    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
		    my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0; 
			    my $ahc = computeMaxProdHitCount($productInfo,$location,
							     $distinctReads,$parameters);
		    my $afh = computeMaxProdFivePrimeHet($productInfo,$parameters);
		    my $pbp = computeProductBasePairing($center,
							$productInfo,$basePairs,$parameters);
		    my $sameShift = computeMaxSameShift($location,
							$distinctReads->{$strand},$productInfo,$parameters);
		    my $bothShift = computeMaxBothShift($basePairs,$location,
							$distinctReads->{$strand},$productInfo,$parameters);
		    # total is the read count for the whole hairpin
		    my $totalSense = 0;
		    foreach my $product (@{$productInfo}) {
			my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			my $length = $relStop - $relStart + 1;
			my $productSequence = substr($sequence,$relStart,$length);
			$totalSense += $adjProdCount;
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
			print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($sample) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$prodList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
			    print HDR "$offset\t$relStart\t$relStop\t$seq";
			    foreach my $bamElement (@{$bamList}) {
				my($sample) = @{$bamElement};
				printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
			    }
			    print HDR "\n";
			}				
		    }
		    #############################
		    # ANTISENSE PRODUCTS
		    my $totalAntisense = 0;
		    foreach my $product (@{$revProductInfo}) {
			my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			$newType = "as-".$newType;
			$totalAntisense += $adjProdCount;
			my $length = $relStop - $relStart + 1;
			my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
			print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($sample) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$prodList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
			    print HDR "$offset\t$relStart\t$relStop\t$seq";
			    foreach my $bamElement (@{$bamList}) {
				my($sample) = @{$bamElement};
				printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
			    }
			    print HDR "\n";
				}				
		    }
		    #############################
		    #############################
		    my($leftCenter, $rightCenter) = @{$center};
		    print HL "$newId\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
		    print HL "$sequence\t$fold\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
			   $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift);
		    
	    
		} else {
		    print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",0,0,$mfe);
		    print HL "$sequence\t$fold\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
			   0,0,0,0,0,0,0,0);		
		}
	    }	   
	}
    }
    close(RRF);
    close(HPL);
    close(HDR);
    $mu->dump();
}

