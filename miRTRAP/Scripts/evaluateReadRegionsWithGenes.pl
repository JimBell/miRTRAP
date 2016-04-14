#!/usr/bin/perl -w
use lib "/home/dhendrix/PerlModules/";
use miRTRAP;
use strict;
$|=1;

my $parameters = {
    "aapdLimit" => 3,
    "ahcLimit" => 5,
    "afhLimit" => 0.5,
    "minSameShift" => 7,
    "minBothShift" => 7,
    "bpDensityLimit" => 0.6,
    "NEIGHBOR" => 1,
    "CHECK_EXONS" => 0,
    "nonMirNeighborLimit" => 10,
    "neighborWindow" => 1000,
    "minGap" => 3,
    "minOverlap" => 4,
    "minMirCount" => 1,
    "minLocusCount" => 5
    };

my $usage = "Usage:\n$0 <config file> <gene models> [valid list] [known miRs]\n";

my $configFile = $ARGV[0] or die $usage;
my $geneModelsFile = $ARGV[1] or die $usage;
my $validList = $ARGV[2] ? $ARGV[2] : "";
my $knownMirsFile = $ARGV[3] ? $ARGV[3] : "";

miRTRAP::readConfigFile($configFile,$parameters);
my $hairpinListFile = $parameters->{hairpinsFile} or die "hairpinListFile: parameter not loaded.\n";
my $hairpinProductsFile = $parameters->{productFile} or die "hairpinProductsfile: parameter not loaded.\n";
my $readRegionsFile = $parameters->{allReadRegionsFile} or die "allReadRegionsFile: parameter not loaded.\n";

my $hairpinList = readHairpinListFile($hairpinListFile);
my $hairpinProducts = readHairpinProductFile($hairpinProductsFile);
my $readRegions = readReadRegions($readRegionsFile);
my $geneModels = readGeneModelFile($geneModelsFile);
my $introns = getIntronsFromGeneModels($geneModels);

if($validList) {
    my $valids = readValidList($validList);
    my $knownMirs = "";
    if($knownMirsFile) {
	$knownMirs = readKnownMirFile($knownMirsFile);		
	# first extract the curatedPositive. This will be used to check against below
	my($curatedPositive,$curatedNegative) = extractCuration($hairpinList,$valids);
	my($curatedPositiveList,$curatedNegativeList) = extractHairpinList($hairpinList,$curatedPositive);
	my($negativeReadRegions) = extractNegativeReadRegions($readRegions,$curatedPositiveList);
	# check the knownMirs against the curatedPositive for contextual data
        printContextualDataWithCuration($knownMirs,$negativeReadRegions,$curatedNegative,"known",$introns,$parameters);
    } else {
	# evaluate the predictedpositive list   
	my($predictedPositive,$predictedNegative,$reasons) = evaluateReadRegions($hairpinList,$hairpinProducts,$geneModels,$introns,$parameters);
	# extract the positive and negative lists
	my($positiveList,$negativeList) = extractHairpinList($hairpinList,$predictedPositive);
	# grab the readRegions that don't overlap positives.
	my($negativeReadRegions) = extractNegativeReadRegions($readRegions,$positiveList);
	my $uniqueHairpinList = extractUniqueHairpinList($positiveList,$hairpinProducts);
	if($parameters->{NEIGHBOR}) {
	    my($finalPositives,$finalNegatives) = collectFinalPredictions($uniqueHairpinList,$negativeReadRegions,$negativeList,$introns,$reasons,$parameters);
	    printCuratedNegatives($finalNegatives,"predicted_vs_curated");
	    evaluatePredictions($finalPositives,$finalNegatives,$valids,$reasons,"predicted_vs_curated");
	} else {
	    printCuratedNegatives($negativeList,"predicted_vs_curated");
	    evaluatePredictions($uniqueHairpinList,$negativeList,$valids,$reasons,"predicted_vs_curated"); 
	}
    }
} else {
    # evaluate the predictedpositive list   
    my($predictedPositive,$predictedNegative,$reasons) = evaluateReadRegions($hairpinList,$hairpinProducts,$geneModels,$introns,$parameters);
    # extract the positive and negative lists
    my($positiveList,$negativeList) = extractHairpinList($hairpinList,$predictedPositive);
    # grab the readRegions that don't overlap positives.
    my($negativeReadRegions) = extractNegativeReadRegions($readRegions,$positiveList);
    my $uniqueHairpinList = extractUniqueHairpinList($positiveList,$hairpinProducts);
    my($finalPositives,$finalNegatives) = collectFinalPredictions($uniqueHairpinList,$negativeReadRegions,$negativeList,$introns,$reasons,$parameters);
    printCuratedNegatives($finalNegatives,"predicted");
    printPredictions($finalPositives,$finalNegatives);
}

sub extractUniqueHairpinList {
    my($hairpinList,$hairpinProducts) = @_;
    my %uniqueHairpinList;
    foreach my $chrom (keys %{$hairpinList}) {
	# sort by total reads, greatest to lowest.
	my @sortedHairpinList = sort {$b->[6] <=> $a->[6]} @{$hairpinList->{$chrom}};
	my @uniqueList;
	foreach my $foldInfo (@sortedHairpinList) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,
	       $aapd,$tapd,$urf,$ahc,$afh,$valid,$shifted) = @{$foldInfo};               
	    my $mirs = getMirs($hairpinProducts->{$tag});
	    unless(insideList($start,$stop,$strand,$mirs,$tag,\@uniqueList)) {
		insertEntry(\@uniqueList,[$start,$stop,$strand,$mirs,$tag,$foldInfo]);	
	    }
	}
	foreach my $uniqueFold (@uniqueList) {
	    my($start,$stop,$strand,$mirs,$tag,$foldInfo) = @{$uniqueFold};
	    push(@{$uniqueHairpinList{$chrom}},$foldInfo);
	}
    }
    return \%uniqueHairpinList
}

sub getMirs {
    my($hairpinProducts) = @_;
    my %mirs;
    foreach my $product (@{$hairpinProducts}) {
	my($side,$type,$prodCount,$maxProdCount,$relStart,$relStop,$seq) = @{$product};	
	if($type eq "miR") {
	    $mirs{$side} = $seq;
	}
    }
    return \%mirs;
}

sub insideList {
    # must overlap at least 50% of max region.
    my($start,$stop,$strand,$mirs,$tag,$list) = @_;
    if(@{$list}) {
	my($searchStart,$searchStop) = getSearchStart($start,$stop,$list);
	# for region in the list:
	for(my $j=$searchStart;$j<=$searchStop;$j++) {
	    my($thisStart,$thisStop,$thisStrand,$theseMirs,$thisTag) = @{$list->[$j]};
	    # check to see if it overlaps the given coords.
	    if(overlap($start,$stop,$strand,$thisStart,$thisStop,$thisStrand)) {				
		if(matchingProducts($mirs,$theseMirs)) {
		    return 1;
		}
	    }
 	}
    }
    return 0;
}

sub getSearchStart {
    my($start,$stop,$list) = @_;
    # searchStart and searchStop are INDICES for the repeat regions in this context
    my $searchStart = 0;
    my $searchStop = scalar(@{$list})-1;    
    return binarySearch($list,$searchStart,$searchStop,$start,$stop);
    #return ($searchStart,$searchStop);
}

sub binarySearch {
    my($list,$searchStart,$searchStop,$start,$stop) = @_;    
    # this 200 is arbitrary. I don't trust it...
    # if the separation is significantly large, try and shorten it.
    if($searchStop - $searchStart > 200) {
        my $searchMiddle = int(($searchStart+$searchStop)/2);
        my($fStart,$fStop) = @{$list->[$searchStart]};
        my($lStart,$lStop) = @{$list->[$searchStop]};
        my($mStart,$mStop) = @{$list->[$searchMiddle]};
        if(($fStart <= $start)&&($stop <= $mStop)) {
            return binarySearch($list,$searchStart,$searchMiddle,$start,$stop);
        } elsif(($mStart  <= $start)&&($stop <= $lStop)) {
            return binarySearch($list,$searchMiddle,$searchStop,$start,$stop);
        } else {
            return ($searchStart,$searchStop);
        }
    } else {
        return ($searchStart,$searchStop);
    }
}

sub matchingProducts {
    my($mirs,$theseMirs) = @_;
    my $MATCH = 0;
    foreach my $side (keys %{$theseMirs}) {
	if($mirs->{$side}) {
	    if($mirs->{$side} eq $theseMirs->{$side}) {
		$MATCH = 1;;
	    } else {
		return 0;
	    }
	}
    }
    return $MATCH;
}

sub overlap {
    my($start,$stop,$strand,$thisStart,$thisStop,$thisStrand) = @_;
    if($strand eq $thisStrand) {
	if(lociContained($start,$stop,$thisStart,$thisStop)) {
	    return 1;
	}
	my $overlap = getOverlap($start,$stop,$thisStart,$thisStop);
	if($overlap) {
	    my $maxLen = max($stop-$start,$thisStop-$thisStart);
	    if($overlap/$maxLen > 0.5) {
		return 1;
	    }
	}
    }
    return 0;
}

sub getOverlap {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    if(($relStart1 <= $relStart2)&&($relStart2 <= $relStop1)) {
	# relStart2 within first window
	return $relStop1 - $relStart2 + 1;
    } elsif(($relStart1 <= $relStop2)&&($relStop2 <= $relStop1)) {
	# if relStop2 within first window
	return $relStop2 - $relStart1 + 1;	
    } elsif(($relStart2 <= $relStart1)&&($relStart1 <= $relStop2)) {
	return $relStop2 - $relStart1 + 1;
    } elsif(($relStart2 <= $relStop1)&&($relStop1<=$relStop2)) {
	return $relStop2-$relStop1 + 1;
    }
    return 0;
}

sub lociContained {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    if(($relStart1 <= $relStart2)&&($relStop2 <= $relStop1)) {
	return 1;
    } elsif(($relStart2 <= $relStart1)&&($relStop1 <= $relStop2)) {
	return 1;
    }
}

sub insertEntry {
    my($list,$entry) = @_;
    if(@{$list}) {
	my $i = @{$list}-1;
	while(($i >= 0)&&($list->[$i]->[0] > $entry->[0])) {
	    $list->[$i+1] = $list->[$i];
	    $i--;
	}
	$list->[$i+1] = $entry;	
    } else {
	push(@{$list},$entry);
    }
}


sub max {
    my @array = @_;
    my @sortedArray = sort {$a <=> $b} @array;
    my $max = pop(@sortedArray);
    return $max;
}

sub min {
    my @array = @_;
    my @sortedArray = sort {$a <=> $b} @array;
    my $min = shift(@sortedArray);
    return $min;
}

sub evaluatePredictions {
    my($positiveList,$negativeList,$valids,$reasons,$label) = @_;
    open(TP,">".$label."_truePos.txt");
    open(FP,">".$label."_falsePos.txt");
    open(FN,">".$label."_falseNeg.txt");
    open(FNTL,">".$label."_falseNeg.tagList");    
    open(PPL,">predicted_positives.txt");
    foreach my $chrom (keys %{$positiveList}) {
	foreach my $hairpin (@{$positiveList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$leftCenter,$rightCenter,$total,$mfe,$seq,$fold,
	       $aapd,$tapd,$urf,$ahc,$afh,$bpDensity,$sameShift,$bothShift) = @{$hairpin};
	    my $location = "$chrom:$start..$stop";
	    if(overlapsList($chrom,$start,$stop,$valids)) {
		# true positive prediction
		print TP "$tag\t$location\t$strand\t$total\n";
	    } else {
		# false positive prediction
		print FP "$tag\t$location\t$strand\t$total\n";
	    }
	    print PPL "$tag\t$chrom:$start..$stop:$strand\t$leftCenter\t$rightCenter\t$mfe\t$seq\t$fold\n";
	}
    }
    my @positiveSearch;
    my @negativeSearch;	
    foreach my $chrom (keys %{$positiveList}) {
	foreach my $hairpin (@{$positiveList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$leftCenter,$rightCenter) = @{$hairpin};
	    push(@positiveSearch,[$chrom,$start,$stop,$strand,$tag]);
	}
    }
    foreach my $chrom (keys %{$negativeList}) {
	foreach my $hairpin (@{$negativeList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$leftCenter,$rightCenter) = @{$hairpin};
	    push(@negativeSearch,[$chrom,$start,$stop,$strand,$tag]);
	}
    }    
    foreach my $term (@{$valids}) {
	my($chrom,$start,$stop,$strand,$tag,$totalReads) = @{$term};
	my $location = "$chrom:$start..$stop";
	unless(overlapsList($chrom,$start,$stop,\@positiveSearch)) {
	    # false negative prediction
	    my $length = $stop - $start + 1;
	    if(my @idList = overlapsList($chrom,$start,$stop,\@negativeSearch)) {
		my $reason = "";;
		foreach my $id (@idList) {
		    $reason .= "$id\t".$reasons->{$id} . "\t";
		    print FNTL "$id\n";
		}
		# false negative prediction		
		print FN "$tag\t$location\t$strand\t$totalReads\t$reason\n";
	    } else {
		# false negative prediction
		print FN "$tag\t$location\t$strand\t$totalReads\tlost?\n";
	    }
	}
    }
    close(FN);
    close(FNTL);
}

sub printPredictions {
    my($positiveList,$negativeList,$valids,$reasons,$label) = @_;
    open(PPL,">predicted_positives.txt");
    foreach my $chrom (keys %{$positiveList}) {
	foreach my $hairpin (@{$positiveList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$leftCenter,$rightCenter,$total,$mfe,$seq,$fold) = @{$hairpin};
	    print PPL "$tag\t$chrom:$start..$stop:$strand\t$leftCenter\t$rightCenter\t$mfe\t$seq\t$fold\n";
	}
    }
}

sub extractNegativeReadRegions {
    my($readRegions,$curatedPositiveList) = @_;
    my %negativeReadRegions;
    foreach my $chrom (keys %{$readRegions}) {
	foreach my $region (@{$readRegions->{$chrom}}) {
	    my($rId,$rStart,$rStop,$rStrand) = @{$region};
	    my $FOUND=0;
	    if($curatedPositiveList->{$chrom}) {
		foreach my $hp (@{$curatedPositiveList->{$chrom}}) {
		    my($tag,$start,$stop) = @{$hp};
		    if(($start <= $rStart)&&($rStop <= $stop)) {
			$FOUND=1;
		    }
		}
	    }
	    unless($FOUND) {
		push(@{$negativeReadRegions{$chrom}},$region);
	    }
	}
    }
    return \%negativeReadRegions;
}

sub combinePredictions {
    my($positives,$negatives) = @_;
    my %combined;
    foreach my $chrom (keys %{$positives}) {
	foreach my $hp (@{$positives->{$chrom}}) {
	    push(@{$combined{$chrom}},$hp);
	}
    }
    foreach my $chrom (keys %{$negatives}) {
	foreach my $hp (@{$negatives->{$chrom}}) {
	    push(@{$combined{$chrom}},$hp);
	}
    }   
    return \%combined;
}

sub printCuratedNegatives {
    my($hairpinList,$hairpinLabel) = @_;
    open(CN,">$hairpinLabel"."_Negatives.gff");
    foreach my $chrom (keys %{$hairpinList}) {
	for(my $i1=0;$i1<@{$hairpinList->{$chrom}};$i1++) {
	    my($tag1,$start1,$stop1,$strand1) = @{$hairpinList->{$chrom}[$i1]};
	    print CN "$chrom\tcuration\tmatch\t$start1\t$stop1\t.\t$strand1\t.\tID=\"$tag1\";\n";
	}
    }
}

sub printCuratedValues {
    my($positiveList,$negativeList,$hairpinLabel) = @_;
    open(PV,">$hairpinLabel"."_PosValues.txt");
    open(NV,">$hairpinLabel"."_NegValues.txt");
    print PV "tag\ttotal\tmfe\taapd\tahc\tafh\tbp\tsameShift\tbothShift\n";
    foreach my $chrom (keys %{$positiveList}) {
	for(my $i=0;$i<@{$positiveList->{$chrom}};$i++) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,$sequence,$fold,
	       $aapd,$tapd,$urf,$ahc,$afh,$bpDensity,$sameShift,$bothShift) = @{$positiveList->{$chrom}[$i]};
	    print PV "$tag\t$total\t$mfe\t$aapd\t$ahc\t$afh\t$bpDensity\t$sameShift\t$bothShift\n";
	}
    }
    foreach my $chrom (keys %{$negativeList}) {
	for(my $i=0;$i<@{$negativeList->{$chrom}};$i++) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,$sequence,$fold,
	       $aapd,$tapd,$urf,$ahc,$afh,$bpDensity,$sameShift,$bothShift) = @{$negativeList->{$chrom}[$i]};
	    print NV "$tag\t$total\t$mfe\t$aapd\t$ahc\t$afh\t$bpDensity\t$sameShift\t$bothShift\n";
	}
    }
    close(PV);
    close(NV);
}

sub printContextualDataWithCuration {
    my($positiveList,$negativeList,$negativeHash,$hairpinLabel,$introns,$parameters) = @_;
    my $nonMirNeighborLimit = $parameters->{nonMirNeighborLimit} or die "nonMirNeighborLimit: parameter not loaded.(printContextualDataWithCuration)\n";
    my $neighborWindow = $parameters->{neighborWindow} or die "neighborWindow: parameter not loaded.\n";
    open(VV,">$hairpinLabel"."_Values.txt");
    open(FP,">$hairpinLabel"."_falsePos.txt");
    open(TP,">$hairpinLabel"."_truePos.txt");
    open(NTN,">$hairpinLabel"."_trueNeg.txt");
    open(FN,">$hairpinLabel"."_falseNeg.txt");
    foreach my $chrom (keys %{$positiveList}) {
	for(my $i1=0;$i1<@{$positiveList->{$chrom}};$i1++) {
	    my($tag1,$start1,$stop1,$strand1,$hStart1,$hStop1) = @{$positiveList->{$chrom}[$i1]};
	    my($leftBound,$rightBound) = getIntronBounds($introns,$chrom,$positiveList->{$chrom}[$i1],$parameters);
	    #print "$tag1 $chrom $start1 $stop1 $strand1\n";
	    my $negativeCount = 0;
	    my $location = "$chrom:$start1..$stop1";
	    if($negativeList->{$chrom}) {
		for(my $i2=0;$i2<@{$negativeList->{$chrom}};$i2++) {
		    my($tag2,$start2,$stop2,$strand2) = @{$negativeList->{$chrom}[$i2]};
		    if(($leftBound <= $start2)&&($stop2 <= $rightBound)) {
			$negativeCount++;
		    }
		}
	    }
	    my $thisNonMirNeighborLimit = int($nonMirNeighborLimit*($rightBound-$leftBound+1)/($stop1-$start1+1+2*$neighborWindow)+0.4999999999);
	    $thisNonMirNeighborLimit = $thisNonMirNeighborLimit > 0 ? $thisNonMirNeighborLimit : 1;
	    print VV "$tag1\t$location\t$strand1\t$negativeCount\t$thisNonMirNeighborLimit ($leftBound,$rightBound)\n";
	    if($negativeCount <= $thisNonMirNeighborLimit) {
		# these examples should be all predicted positive.
		if($negativeHash->{$tag1}) {
		    # this is in the negative curated set, therefore false positive
		    print FP "$tag1\t$location\t$strand1\t$negativeCount\n";
		} else {
		    # this is positive in the curated set. therefore true positive.
		    print TP "$tag1\t$location\t$strand1\t$negativeCount\n";
		}
	    } else {
		# these examples should be negative, since they have too many negative neighbors.
		if($negativeHash->{$tag1}) {
                    # this is in the negative curated set, therefore true negative by neighbor rule.
		    print NTN "$tag1\t$location\t$strand1\t$negativeCount\n";
                } else {
                    # this is positive in the curated set, but now predicted negative
                    print FN "$tag1\t$location\t$strand1\t$negativeCount\tFails non-miR neighbor count\n";
                }
	    }
	}
    }
}

sub collectFinalPredictions {
    my($positiveList,$negativeRegions,$negativeList,$introns,$reasons,$parameters) = @_;
    my $nonMirNeighborLimit = $parameters->{nonMirNeighborLimit} or die "nonMirNeighborLimit: parameter not loaded. (collectFinalPredictions)\n";
    my $neighborWindow = $parameters->{neighborWindow} or die "neighborWindow: parameter not loaded.\n";
    my %finalPrediction;
    foreach my $chrom (keys %{$positiveList}) {
	for(my $i1=0;$i1<@{$positiveList->{$chrom}};$i1++) {
	    my($tag1,$start1,$stop1,$strand1,$hStart1,$hStop1,$total1,$mfe1,$seq1,$fold1) = @{$positiveList->{$chrom}[$i1]};
	    my($leftBound,$rightBound) = getIntronBounds($introns,$chrom,$positiveList->{$chrom}[$i1],$parameters);
	    #print "$tag1 $chrom $start1 $stop1 $strand1\n";
	    my $negativeCount = 0;
	    my $location = "$chrom:$start1..$stop1";
	    if($negativeRegions->{$chrom}) {
		for(my $i2=0;$i2<@{$negativeRegions->{$chrom}};$i2++) {
		    my($tag2,$start2,$stop2,$strand2) = @{$negativeRegions->{$chrom}[$i2]};
		    if(($leftBound <= $start2)&&($stop2 <= $rightBound)) {
			$negativeCount++;
		    }
		}
	    }
	    my $thisNonMirNeighborLimit = int($nonMirNeighborLimit*($rightBound-$leftBound+1)/($stop1-$start1+1+2*$neighborWindow)+0.4999999999);
	    $thisNonMirNeighborLimit = $thisNonMirNeighborLimit > 0 ? $thisNonMirNeighborLimit : 1;
	    if($negativeCount <= $thisNonMirNeighborLimit) {
	    # these examples should be all predicted positive.		
		push(@{$finalPrediction{$chrom}},$positiveList->{$chrom}[$i1]);
	    } else {
		print "$tag1\tFailed non-miR neighbor count test. $negativeCount vs $thisNonMirNeighborLimit\n";
		push(@{$negativeList->{$chrom}},$positiveList->{$chrom}[$i1]);
		$reasons->{$tag1} = "Failed non-miR neighbor count test. $negativeCount vs $thisNonMirNeighborLimit.";
	    }
	}
    }
    return(\%finalPrediction,$negativeList);
}

sub appendFalseNegatives {
    my($lostCuration,$negativeList,$negativeHash,$hairpinLabel,$introns,$parameters) = @_;
    my $nonMirNeighborLimit = $parameters->{nonMirNeighborLimit} or die "nonMirNeighborLimit: parameter not loaded.(appendFalseNegatives)\n";
    my $neighborWindow = $parameters->{neighborWindow} or die "neighborWindow: parameter not loaded.\n";
    open(VV,">>$hairpinLabel"."_Values.txt");
    open(FN,">>$hairpinLabel"."_falseNeg.txt");
    foreach my $chrom (keys %{$lostCuration}) {
	for(my $i1=0;$i1<@{$lostCuration->{$chrom}};$i1++) {
	    my($tag1,$start1,$stop1,$strand1,$hStart,$hStop,$reason) = @{$lostCuration->{$chrom}[$i1]};	    
	    my($leftBound,$rightBound) = getIntronBounds($introns,$chrom,$lostCuration->{$chrom}[$i1],$parameters);	    
	    my $negativeCount = 0;
	    my $location = "$chrom:$start1..$stop1";
	    if($negativeList->{$chrom}) {
		for(my $i2=0;$i2<@{$negativeList->{$chrom}};$i2++) {
		    my($tag2,$start2,$stop2,$strand2) = @{$negativeList->{$chrom}[$i2]};
		    if(($leftBound <= $start2)&&($stop2 <= $rightBound)) {
			$negativeCount++;
		    }
		}
	    }
	    my $thisNonMirNeighborLimit = int($nonMirNeighborLimit*($rightBound-$leftBound+1)/($stop1-$start1+1+2*$neighborWindow)+0.4999999999);
	    $thisNonMirNeighborLimit = $thisNonMirNeighborLimit > 0 ? $thisNonMirNeighborLimit : 1;
	    print VV "$tag1\t$location\t$strand1\t$negativeCount\t$thisNonMirNeighborLimit ($leftBound,$rightBound)\n";
	    print FN "$tag1\t$location\t$strand1\t$negativeCount\t$reason\n";
	}
    }
}

sub getIntronBounds {
    my($introns,$chrom,$hairpin,$parameters) = @_;
    my $neighborWindow = $parameters->{neighborWindow} or die "neighborWindow: parameter not loaded.\n";
    my($tag,$start,$stop,$strand,$hStart,$hStop) = @{$hairpin};    
    my @intronList;
    my $thisStart = $start + 50;
    my $thisStop = $stop - 50;
    my $leftBound = $thisStart - $neighborWindow;
    my $rightBound = $thisStop + $neighborWindow;
    $leftBound = $leftBound > 0 ? $leftBound : 0;
    #print "starting with $leftBound $rightBound\n";
    if($introns->{$chrom}{$strand}) {	
	foreach my $intron (@{$introns->{$chrom}{$strand}}) {
	    my($inStart,$inStop,$geneId,$type1,$type2) = @{$intron};
	    if(($inStart<=$thisStart+5)&&($stop-5<=$inStop)) {
		push(@intronList,[$inStart,$inStop,$geneId]);
	    }
	}
    }
    if(@intronList) {
	if(@intronList == 1) {
	    my $intron = shift(@intronList);
	    my($inStart,$inStop) = @{$intron};
	    #print "found: $inStart $inStop\n";
	    if($thisStart - $inStart <= $neighborWindow) {
		$leftBound = $inStart;
	    }
	    if($inStop - $stop <= $neighborWindow) {
		$rightBound = $inStop;
	    }
	    #print "loaded: $leftBound $rightBound\n";
	    return ($leftBound,$rightBound);
	} else {
	    my $minLength = $rightBound - $leftBound;
	    foreach my $intron (@intronList) {
		my($inStart,$inStop,$geneId) = @{$intron};
		if($inStop - $inStart < $minLength) {
		    $minLength = $inStop - $inStart;
		    $leftBound = $inStart;
		    $rightBound = $inStop;
		}		
	    }
	    return ($leftBound,$rightBound);
	}
    } else {
	return ($leftBound,$rightBound);
    }
}

sub overlapsIntron {
    my($chrom,$hairpin,$introns) = @_;
    my($tag,$start,$stop,$strand,$hStart,$hStop) = @{$hairpin};
    my $thisStart = $start - 50;
    my $thisStop = $stop - 50;
    if($introns->{$chrom}{$strand}) {	
	if($introns->{$chrom}{$strand}) {
	    foreach my $intron (@{$introns->{$chrom}{$strand}}) {
		my($inStart,$inStop,$geneId,$type1,$type2) = @{$intron};
		if(($inStart<=$thisStart+5)&&($thisStop-5<=$inStop)) {
		    return 1;
	    }
	    }
	}
    }
    return 0;
}

sub overlapsExon {
    my($chrom,$hairpin,$exons) = @_;
    my($tag,$start,$stop,$strand,$hStart,$hStop) = @{$hairpin};
    my $thisStart = $start + 50;
    my $thisStop = $stop - 50;
    foreach my $geneId (keys %{$geneModels}) {
	foreach my $exon (@{$geneModels->{$geneId}}) {
	    my($type1,$chrom1,$start1,$stop1,$strand1) = @{$exon};	    
	    if($chrom eq $chrom1) {
		if($strand eq $strand1) {
		    if(($start1 <= $thisStart)&&($thisStop <= $stop1)) {
			return 1;
		    }
		}
	    }
	}
    }
    return 0;
}

sub extractHairpinList {
    my($hairpinList,$predictedPositive) = @_;
    my %positiveList;
    my %negativeList;
    foreach my $chrom (keys %{$hairpinList}) {
	foreach my $hairpin (@{$hairpinList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,$sequence,$fold,
	       $aapd,$tapd,$urf,$ahc,$afh,$bpDensity,$sameShift,$bothShift) = @{$hairpin};
	    if($predictedPositive->{$tag}) {
		push(@{$positiveList{$chrom}},$hairpin);
	    } else {
		push(@{$negativeList{$chrom}},$hairpin);
	    }
	}
    }
    return(\%positiveList,\%negativeList);
}

sub extractCuration {
    my($hairpinList,$valids) = @_;
    my %curatedPositive;
    my %curatedNegative;
    foreach my $chrom (keys %{$hairpinList}) {
	foreach my $hairpin (@{$hairpinList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop) = @{$hairpin};
	    my $thisStart = $start + 50;
	    my $thisStop = $stop - 50;
	    if(overlapsList($chrom,$thisStart,$thisStop,$valids)) {
		$curatedPositive{$tag}++;
	    } else {
		$curatedNegative{$tag}++;
	    }
	}
    }
    return (\%curatedPositive,\%curatedNegative);
}

sub extractLostCuration {
    my($positiveList,$valids,$negativeList,$reasons) = @_;
    my @positiveSearch;
    my @negativeSearch;
    foreach my $chrom (keys %{$positiveList}) {
	foreach my $hairpin (@{$positiveList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop) = @{$hairpin};
	    push(@positiveSearch,[$chrom,$start,$stop,$strand,$tag]);
	}
    }
    foreach my $chrom (keys %{$negativeList}) {
	foreach my $hairpin (@{$negativeList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop) = @{$hairpin};
	    push(@negativeSearch,[$chrom,$start,$stop,$strand,$tag]);
	}
    }    
    my %lostCuration;
    foreach my $term (@{$valids}) {
	my($chrom,$start,$stop,$strand,$id) = @{$term};
	#print "checking $id $chrom $start $stop\n";
	unless(overlapsList($chrom,$start,$stop,\@positiveSearch)) {
	    my $length = $stop - $start + 1;
	    if(my @idList = overlapsList($chrom,$start,$stop,\@negativeSearch)) {
		my $reason = "";;
		foreach my $id (@idList) {
		    $reason .= "$id\t".$reasons->{$id} . "\t";
		}
		push(@{$lostCuration{$chrom}},[$id,$start,$stop,$strand,0,$length,$reason]);
	    } else {
		push(@{$lostCuration{$chrom}},[$id,$start,$stop,$strand,0,$length,"lost"]);	
	    }
	}
    }
    return \%lostCuration;
}

sub evaluateReadRegions {
    my($hairpinList,$hairpinProducts,$geneModels,$introns,$parameters) = @_;
    my %predictedPositive;
    my %predictedNegative;
    my %reasons;
    foreach my $chrom (keys %{$hairpinList}) {
	foreach my $hairpin (@{$hairpinList->{$chrom}}) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,$sequence,$fold,
	       $aapd,$tapd,$urf,$ahc,$afh,$bpDensity,$sameShift,$bothShift) = @{$hairpin};	    
	    my($accepted,$reason) = evaluate($chrom,$hairpin,$hairpinProducts,$parameters);
	    if($accepted) {
		$predictedPositive{$tag}++;
	    } else {
		$predictedNegative{$tag}++;
		print "$tag\t$reason\n";
		$reasons{$tag} = $reason;
	    }
	}
    }
    return (\%predictedPositive,\%predictedNegative,\%reasons);
}

sub evaluate {
    my($chrom,$hairpin,$hairpinProducts,$parameters) = @_;
    my $minSameShift = $parameters->{minSameShift} or die "minSameShift: parameter not loaded.\n";
    my $minBothShift = $parameters->{minBothShift} or die "minBothShift: parameter not loaded.\n";
    my $bpDensityLimit = $parameters->{bpDensityLimit} or die "bpDensityLimit: parameter not loaded.\n";
    my $aapdLimit = $parameters->{aapdLimit} or die "aapdLimit: parameter not loaded.\n";
    my $ahcLimit = $parameters->{ahcLimit} or die "ahcLimit: parameter not loaded.\n";
    my $afhLimit = $parameters->{afhLimit} or die "afhLimit: parameter not loaded.\n";
    my $minMirCount = $parameters->{minMirCount} or die "minMirCount: parameter not loaded.\n";
    my $minLocusCount = $parameters->{minLocusCount} or die "minLocusCount: parameter not loaded.\n";
    
    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,$sequence,$fold,
       $aapd,$tapd,$urf,$ahc,$afh,$bpDensity,$sameShift,$bothShift) = @{$hairpin};
    my $reason;
    #my $asymmetry = computeAsymmetryScore($hairpin);
    if($sameShift > $minSameShift) {
	$reason .= "Has shifted products on the same arm.";
    }
    if($bothShift > $minBothShift) {
	$reason .= "Has shifted products on opposite arms.";
    }    
    if($bpDensity < $bpDensityLimit) {
	$reason .= "miR outside of bp limit: $bpDensity.";
    }	
    if($aapd > $aapdLimit) {
	$reason .= "Fails the antisense displacement test: $aapd. ";
    }
    if($ahc > $ahcLimit) {
	$reason .= "Fails average hit count test: $ahc. ";
    }
    if($afh >= $afhLimit) {
	$reason .= "Fails 5' Hetergeneity test: $afh. ";
    }
    if($total < $minLocusCount) {
	$reason .= "Fails because total reads = $total < $minLocusCount = minLocusCount. ";
    }
    if($hairpinProducts->{$tag}) {
	if(hasSplitProducts($hairpinProducts->{$tag},$parameters)) {
	    $reason .= "Has products spanning loop and arm. ";
	}
	if(hasSpacedProducts($tag,$total,$hairpinProducts,$parameters)) {
	    $reason .= "Has discontinuous products, unlikely dicer cut. ";
	}
	if(hasOverlappingProducts($tag,$total,$hairpinProducts,$parameters)) {
	    $reason .= "Has overlapping products, unlikely dicer cut. ";
	}
	if(lowAbundantMirs($hairpinProducts->{$tag},$parameters)) {
	    $reason .= "Does not have a miR > $minMirCount read. ";
	}
    } else {
	$reason .= "Does not have products. ";
    }
    unless($reason) {
	if($parameters->{CHECK_EXONS}) {
	    unless(overlapsIntron($chrom,$hairpin,$introns,$parameters)) {
		if(overlapsExon($chrom,$hairpin,$geneModels,$parameters)) {
		    $reason .= "Overlaps exon. ";
		}
	    }
	}
    }
    if($reason) {
	return (0,$reason);
    } else {
	return (1,$reason);
    }
}

sub hasSplitProducts {
    my($products,$parameters) = @_;
    foreach my $product (@{$products}) {
	my($side,$type,$prodCount,$maxProdCount,$relStart,$relStop) = @{$product};
	if($type eq "split") {
	    return 1;
	}
    }
    return 0;
}

sub hasSpacedProducts {
    my($tag,$total,$hairpinProducts,$parameters) = @_;
    my $minGap = $parameters->{minGap} or die "minGap: parameter not loaded.\n";
    my @sortedProducts = sort {$a->[4] <=> $b->[4]} @{$hairpinProducts->{$tag}};
    #print "checking $tag\n";
    for(my $i1=0;$i1<@sortedProducts;$i1++) {
	my $product1 = $sortedProducts[$i1];
	my($side1,$type1,$prodCount1,$maxProdCount1,$relStart1,$relStop1) = @{$product1};
	if($type1 eq "miR") {
	    #print "evaluating $side1 $type1\n";
	    # if there is an 5' product
	    if((0 <= $i1-1)&&($side1 eq "5p")) {
		my $product2 = $sortedProducts[$i1-1];
		my($side2,$type2,$prodCount2,$maxProdCount2,$relStart2,$relStop2) = @{$product2};
		if(($type2 eq fivePrimeNeighborType($side1))&&($prodCount2 >= 1)&&($prodCount2 >= 0.05*$total)) {
		    if($relStart1 - $relStop2 - 1 > $minGap) {
			return 1;
		    }
		}
	    }
	    # if there is a 3' product
	    if(($i1+1 <= $#sortedProducts)&&($side1 eq "3p")) {
		my $product2 = $sortedProducts[$i1+1];
		my($side2,$type2,$prodCount2,$maxProdCount2,$relStart2,$relStop2) = @{$product2};
		if(($type2 eq threePrimeNeighborType($side1))&&($prodCount2 >= 1)&&($prodCount2 >= 0.05*$total)) {
		    #print "comparing to $side2 $type2 ($relStart1,$relStop1) ($relStart2,$relStop2)\n";
		    if($relStart2 - $relStop1 - 1 > $minGap) {
			return 1;
		    }
		}
	    }
	}
    }
    return 0;
}

sub hasOverlappingProducts {
    my($tag,$total,$hairpinProducts) = @_;
    my $minOverlap = $parameters->{minOverlap} or die "minOverlap: parameter not loaded.\n";
    my @sortedProducts = sort {$a->[4] <=> $b->[4]} @{$hairpinProducts->{$tag}};
    #print "checking $tag\n";
    for(my $i1=0;$i1<@sortedProducts;$i1++) {
	my $product1 = $sortedProducts[$i1];
	my($side1,$type1,$prodCount1,$maxProdCount1,$relStart1,$relStop1) = @{$product1};
	if($type1 eq "miR") {
	    #print "evaluating $side1 $type1\n";
	    # if there is an 5' product
	    if((0 <= $i1-1)&&($side1 eq "5p")) {
		my $product2 = $sortedProducts[$i1-1];
		my($side2,$type2,$prodCount2,$maxProdCount2,$relStart2,$relStop2) = @{$product2};
		if(($type2 eq fivePrimeNeighborType($side1))&&($prodCount2 >= 1)&&($prodCount2 >= 0.05*$total)) {
		    if($relStop2 - $relStart1 + 1 >= $minOverlap) {
			return 1;
		    }
		}
	    }
	    # if there is a 3' product
	    if(($i1+1 <= $#sortedProducts)&&($side1 eq "3p")) {
		my $product2 = $sortedProducts[$i1+1];
		my($side2,$type2,$prodCount2,$maxProdCount2,$relStart2,$relStop2) = @{$product2};
		if(($type2 eq threePrimeNeighborType($side1))&&($prodCount2 >= 1)&&($prodCount2 >= 0.05*$total)) {
		    #print "comparing to $side2 $type2 ($relStart1,$relStop1) ($relStart2,$relStop2)\n";
		    if($relStop1 - $relStart2 + 1 >= $minOverlap) {
			return 1;
		    }
		}
	    }
	}
    }
    return 0;
}

sub lowAbundantMirs {
    my($products,$parameters) = @_;
    my $minMirCount = $parameters->{minMirCount} or die "minMirCount: parameter not loaded.\n";
    my $FOUND = 1;
    foreach my $product (@{$products}) {
	my($side,$type,$prodCount,$maxProdCount,$relStart,$relStop) = @{$product};
	if($type eq "miR") {
	    if($prodCount > $minMirCount) {
		$FOUND = 0;
	    }
	}
    }
    return $FOUND;
}

sub computeAsymmetryScore {
    my($hairpin) = @_;
    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,$sequence,$fold,
       $aapd,$tapd,$urf,$ahc,$afh,$bpDensity,$sameShift,$bothShift) = @{$hairpin};
    my(@centers) = miRTRAP::getHairpinCenters($fold);
    foreach my $center (@centers) {
	my($leftCenter,$rightCenter) = @{$center};		    
	if(($hStart <= $leftCenter)&&($rightCenter <= $hStop)) {
	    my $fivePrimeNt = $leftCenter - $hStart + 1;
	    my $threePrimeNt = $hStop - $rightCenter + 1;
	    my $score = abs($fivePrimeNt - $threePrimeNt)/($fivePrimeNt + $threePrimeNt);
	    return $score;
	}
    }
}

sub fivePrimeNeighborType {
    my($side) = @_;
    if($side eq "5p") {
	return "moR";
    } else {
	die "unexpected side in fivePrimeNeighborType: $side\n";
    }
}

sub threePrimeNeighborType {
    my($side) = @_;
    if($side eq "3p") {
	return "moR";
    } else {
	die "unexpected side in threePrimeNeighborType: $side\n";
    }
}


sub overlapsList {
    my($chrom,$start,$stop,$list) = @_;
    my @IDLIST;
    foreach my $item (@{$list}) {
	my($tChrom,$tStart,$tStop,$tStrand,$id) = @{$item};
	if($chrom eq $tChrom) {
	    if((($start <= $tStart)&&($tStart <= $stop))||
	       (($start <= $tStop)&&($tStop <= $stop))||
	       (($tStart <= $start)&&($start <= $tStop))||
	       (($tStart <= $stop)&&($stop <= $tStop))) {
		push(@IDLIST,$id);
	    }
	}
    }
    return @IDLIST;
}

sub readValidList {
    my $listFile = shift;
    open(LF,$listFile) or die "could not open $listFile\n";
    my @locusList;
    while(<LF>) {
	chomp;
	unless(/\#/) {
	    my($id,$location,$strand,$totalReads) = split(/\t/);
	    my($chrom,$start,$stop) = miRTRAP::parseLocation($location);	    
	    push(@locusList,[$chrom,$start,$stop,$strand,$id,$totalReads]);
	}
    }
    return \@locusList;
}

sub readKnownMirFile {
    my $miRNAFile = shift;
    my %knownSet;
    open(MRNAFILE,$miRNAFile);
    while(<MRNAFILE>) {
        unless(/\#/) {
            chomp;
            my($chrom,$source,$type,$start,$stop,$score,$strand,$dot,$id)=split(/\t/,$_);
	    my $uniqueId;
	    my $thisId;
            if(($thisId) = $id =~ /ID=\"(.*)\";/) {
		$uniqueId = $thisId;
	    } elsif(($thisId) = $id =~ /ID=(.*?);/) {
		$uniqueId = $thisId;
	    }
	    my $length = $stop - $start + 1;
	    push(@{$knownSet{$chrom}},[$uniqueId,$start,$stop,$strand,0,$length]);
        }
    }
    return \%knownSet;
}

sub readReadRegions {
    my $readRegionsFile = shift;
    my %readRegions;
    open(RRF,$readRegionsFile) or die "could not open $readRegionsFile\n";
    while(<RRF>) {
	my($id,$location,$length) = split();
	my($chrom,$start,$stop,$strand) = miRTRAP::parseLocation($location);
	push(@{$readRegions{$chrom}},[$id,$start,$stop,$strand]);
    }    
    return \%readRegions;
}

sub readHairpinListFile {
    my $hairpinFile = shift;
    my %hairpinList;
    open(MFILE,$hairpinFile) or die "could not open $hairpinFile\n";
    while(<MFILE>) {
        chomp;
        unless(/\#/) {
	    my($tag,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$total,$mfe,$seq,$fold,
	       $arpd,$trpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift) = split(/\t/);
	    push(@{$hairpinList{$chrom}},[$tag,$start,$stop,$strand,$leftCenter,$rightCenter,$total,$mfe,
					  $seq,$fold,$arpd,$trpd,$urf,$ahc,$afh,$pbp,
					  $sameShift,$bothShift]);
	}
    }
    return \%hairpinList;
}

sub readHairpinProductFile {
    my $hairpinProductFile = shift;
    my %hairpinProducts;
    open(HPF,$hairpinProductFile) or die "could not open $hairpinProductFile\n";
    while(<HPF>) {
	chomp;
	unless(/\#/) {
	    my($tag,$side,$type,$prodCount,$maxProdCount,$relStart,$relStop,$productSequence,@sampleCounts) = split(/\t/);
	    push(@{$hairpinProducts{$tag}},[$side,$type,$prodCount,$maxProdCount,$relStart,$relStop,$productSequence]);
	}
    }
    return \%hairpinProducts;
}

sub readGeneModelFile {
    my($geneModelGffFile) = @_;
    my %geneModels;
    my %mRNA;
    open(FILEA,$geneModelGffFile) or die "could not open $geneModelGffFile";
    while(<FILEA>) {
        unless(/\#/) {
	    chomp;
            my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info)=split(/\t/);
	    if($type eq "mRNA") {
		if($info =~ /ID=(.*?);/) {
		    my($transcriptId) = $info =~ /ID=(.*?);/;
		    $mRNA{$transcriptId}++;
		} elsif($info =~ /ID=(.*)/) {
		    my($transcriptId) = $info =~ /ID=(.*)/;
		    $mRNA{$transcriptId}++;
		}
	    }
	}
    }
    open(FILEA,$geneModelGffFile) or die "could not open $geneModelGffFile";
    while(<FILEA>) {
        unless(/\#/) {
            chomp;
            my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info)=split(/\t/);
	    if($type eq "exon") {
                my($parent) = $info =~ /Parent=(.*);?/;
		my @parents = split(/\,/,$parent);		
		foreach my $geneId (@parents) {		    
		    if($mRNA{$geneId}) {
			push(@{$geneModels{$geneId}},[$type,$chrom,$start,$stop,$strand]);
		    }
		}
	    }
	}
    }
    return \%geneModels;
}

sub getIntronsFromGeneModels {
    my($geneModels) = @_;
    my %introns;
    my %USED;
    foreach my $geneId (keys %{$geneModels}) {
	my @exons;
	foreach my $term (@{$geneModels->{$geneId}}) {
	    my($type,$chrom,$start,$stop,$strand) = @{$term};
	    if($type eq "exon") {
		push(@exons,$term);
	    }
	}
	my @sorted = sort {$a->[2] <=> $b->[2]} @exons;
	for(my $i=0;$i<@sorted-1;$i++) {
	    my($type1,$chrom1,$start1,$stop1,$strand1) = @{$exons[$i]};
	    my($type2,$chrom2,$start2,$stop2,$strand2) = @{$exons[$i+1]};
	    if($stop1+1 < $start2) {
		unless($USED{$chrom1.$stop1.$start2.$strand1}) {
		    push(@{$introns{$chrom1}{$strand1}},[$stop1+1,$start2-1,$geneId,$type1,$type2]);
		    $USED{$chrom1.$stop1.$start2.$strand1}++;
		}
	    }
	}
    }
    return \%introns;
}

