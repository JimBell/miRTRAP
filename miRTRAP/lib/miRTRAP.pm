package miRTRAP;
use strict;
use warnings;

sub readConfigFile {
    my($configFile,$parameters) = @_;
    open(CONFIGFILE,$configFile) or die "FAIL: could not open $configFile\n";
    while(<CONFIGFILE>) {
	chomp;
	my($key,$value) = split(/\s+=\s+/);
	$parameters->{$key} = $value;
    }
    my $filePrefix = $parameters->{filePrefix} or die "FAIL: filePrefix not loaded.\n";
    $parameters->{distancesFile} = $filePrefix."_distances.txt";
    $parameters->{readRegionsFile} = $filePrefix.".txt";
    $parameters->{longReadRegionsFile} = $filePrefix."_long.txt";
    $parameters->{allReadRegionsFile} = $filePrefix."_full.txt";
    $parameters->{readRegionsFastaFile} = $filePrefix . ".fasta";
    $parameters->{rnaFoldOutput} = $filePrefix . ".mfe";
    $parameters->{filteredCandidatesFile} = $filePrefix . "_filteredCandidates.txt";
    $parameters->{hairpinsFile} = $filePrefix . "_hairpins.txt";
    $parameters->{distinctReadsFile} = $filePrefix . "_distinctReads.txt";
    $parameters->{productFile} = $filePrefix . "_products.txt";    
    $parameters->{configFile} = $configFile;
    return $parameters;
}

######################################
# SEQUENCE / GENOME TOOLS            #
######################################

sub reverseComplement {
# Returns the reverse complement of the input sequence.
    my($seq)=@_;
    $seq =~ tr/acgtACGT/tgcaTGCA/;
    $seq=reverse($seq);
    return $seq;
}

sub parseLocation {
    my($location)=@_;
    my($chrom,$start,$end,$strand);
    
    if($location =~ /(.*)\:(-?\d+)\-(-?\d+)\:(.*)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;
    } elsif($location =~ /(.*)\:(-?\d+)\-(-?\d+)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)\:(.*)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;   	
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";   	
    }
    return ($chrom,$start,$end,$strand);
}    

sub loadGenome {
    my $genomeFile = shift;
    my $chromRegExp = shift;
    
    my($chromRegExp1,$chromRegExp2);
    if($chromRegExp) {
	($chromRegExp1,$chromRegExp2) = $chromRegExp =~ /s\/(.*)\/(.*)\// or
	    die "regular expression must be of the form \"s/a/b/\"\nto change from 'a' to 'b'\n";;
    }
    my %sequences;
    my $tag;
    open(GFILE,$genomeFile) or die "could not open $genomeFile\n";
    while(<GFILE>) {
        chomp;
        if(/>/) {
            s/>//g;
            my @terms = split;
	    $tag = shift(@terms);
	    if($chromRegExp) {		
		$tag =~ s/$chromRegExp1/$chromRegExp2/;
	    }
        } else {
            $sequences{$tag} .= $_;
        }
    }
    return \%sequences;
}

sub getSequence {
    my($location,$sequences)=@_;
    my($chrom,$start,$stop,$strand)=parseLocation($location);
    if($sequences->{$chrom}) {
	#print "trying to get ", $stop-$start+1, "\n";
        my $string =  substr($sequences->{$chrom},$start-1,$stop-$start+1);
        if($strand eq "+") {
            return $string;
        } else {
	    return reverseComplement($string);
        }
    } else {
        die "could not find $chrom in function getSequence\nLoaded from $location\n";
    }
}

############################
# BINARY SEARCH TOOLS      #
############################

sub significantlyInsideList {
    # must overlap at least 50% of max region.
    my($start,$stop,$list) = @_;
    if(@{$list}) {
	my($searchStart,$searchStop) = getSearchStart($start,$stop,$list);
	# for region in the list:
	for(my $j=$searchStart;$j<=$searchStop;$j++) {
	    my($thisStart,$thisStop,$foldInfo) = @{$list->[$j]};
	    # check to see if it overlaps the given coords.
	    if(significantOverlap($start,$stop,$thisStart,$thisStop)) {
		return 1;
	    }
 	}
    }
    return 0;
}

sub significantOverlap {
    my($start,$stop,$thisStart,$thisStop) = @_;
    if(lociContained($start,$stop,$thisStart,$thisStop)) {
	return 1;
    }
    my $overlap = getOverlap($start,$stop,$thisStart,$thisStop);
    if($overlap) {
	my $maxLen = max($stop-$start,$thisStop-$thisStart);
	if($overlap/$maxLen > 0.35) {
	    return 1;
	}
    }
    return 0;
}

sub insideList {
    my($start,$stop,$list) = @_;
    #print "testing $start $stop\n";
    #print "size of list: ", scalar(@{$list}), "\n";
    if($list) {
	if(@{$list}) {
	    my($searchStart,$searchStop) = getSearchStart($start,$stop,$list);
	    # for region in the list:
	    #print "searching from $searchStart to $searchStop\n";
	    for(my $j=$searchStart;$j<=$searchStop;$j++) {
		my($thisStart,$thisStop) = @{$list->[$j]};
		# check to see if it overlaps the given coords.
		if((($thisStart <= $start)&&($start <= $thisStop))||
		   (($thisStart <= $stop)&&($stop <= $thisStop))||
		   (($start <= $thisStart)&&($thisStart <= $stop))||
		   (($start <= $thisStop)&&($thisStop <= $stop))) {		
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

##############
# MATH TOOLS #
##############

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

#########################
# FOLD PROCESSING TOOLS #
#########################

sub recenterHairpins {
    my($foldFile,$folds,$genome,$dataSet,$hitCount,$parameters) = @_;    
    my($fileBase) = $foldFile =~ /(.*).mfe/ or die "could not parse filename: $foldFile. Expecting *.mfe\n";
    my $maxHpl = 65;
    open(RNAF,">".$fileBase.".rc.fasta");
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand (keys %{$folds->{$chrom}}) {
            foreach my $foldInfo (@{$folds->{$chrom}{$strand}}) {
                my($start,$stop,$tag,$mfe,$sequence,$fold) = @{$foldInfo};
		if(multipleHairpins($fold,$parameters)) {
		    my @centers = getMergedHairpinCenters($fold,$parameters);
		    my @labels = ("a","b","c","d","e","f");
		    my $COUNT = 0;
		    foreach my $center (@centers) {
			my($leftCenter,$rightCenter) = @{$center};
			my $newStart = $start + $leftCenter - $maxHpl;
			my $newStop = $start + $rightCenter + $maxHpl;
			if($strand eq "-") {
			    $newStop = $stop - $leftCenter + $maxHpl;
			    $newStart = $stop - $rightCenter - $maxHpl;
			}
			$newStart = max($newStart,1);
			$newStop = min($newStop,length($genome->{$chrom}));
			my($hStart,$hStop) = extendHairpinToFullEnds($fold,$leftCenter,$rightCenter);
			# thisStart/thisStop marks the ends of the hairpins
			my $thisStart = $start+$hStart;
			my $thisStop = $start+$hStop;
			if($strand eq "-") {
			    $thisStart = $stop - $hStop;
			    $thisStop = $stop - $hStart;
			}
			my $readCount = getOverlappingReadCount($chrom,$thisStart,$thisStop,$strand,
								$dataSet,$hitCount);
			#print "($leftCenter,$rightCenter) $hStart $hStop -> $thisStart $thisStop : $readCount\n";
			if($readCount > -1) {
			    my $newLocation = $chrom.":".$newStart."..".$newStop.":".$strand;
			    my $sequence = getSequence($newLocation,$genome);
			    my $newId = $tag . $labels[$COUNT];
			    print RNAF ">$newId\t$newLocation\n$sequence\n";
			    $COUNT++;
			}
		    }
		} else {
		    my $leftCenter = getMaxLeftCenter($fold);
		    my $rightCenter = getMaxRightCenter($fold);  
		    my $newStart = $start + $leftCenter - $maxHpl;
		    my $newStop = $start + $rightCenter + $maxHpl;
		    if($strand eq "-") {
			$newStop = $stop - $leftCenter + $maxHpl;
			$newStart = $stop - $rightCenter - $maxHpl;
		    }
		    $newStart = max($newStart,1);
		    $newStop = min($newStop,length($genome->{$chrom}));                 
		    my $newLocation = $chrom.":".$newStart."..".$newStop.":".$strand;
		    my $sequence = getSequence($newLocation,$genome);
		    print RNAF ">$tag\t$newLocation\n$sequence\n";
		}
	    }
	}
    }
    close(RNAF);
}

sub printUniqueHairpins {
    my $folds = shift;
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand (keys %{$folds->{$chrom}}) {
            my @sortedFoldData = sort {$a->[3] <=> $b->[3]} @{$folds->{$chrom}{$strand}};
            my @uniqueList;
            foreach my $foldInfo (@sortedFoldData) {
                my($start,$stop,$tag,$mfe,$sequence,$fold) = @{$foldInfo};
		my $leftCenter = getMaxLeftCenter($fold);
		my $rightCenter = getMaxRightCenter($fold);
                my($hStart,$hStop) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
                my $absStart = $start + $hStart;
                my $absStop = $start + $hStop;
                if($strand eq "-") {
                    $absStart = $stop - $hStop;
                    $absStop = $stop - $hStart;
                }
                unless(significantlyInsideList($absStart,$absStop,\@uniqueList)) {
                    insertEntry(\@uniqueList,[$absStart,$absStop,
					      $start,$stop,$tag,$mfe,$sequence,$fold]);
                }
            }
            foreach my $foldInfo (@uniqueList) {
                my($absStart,$absStop,
                   $start,$stop,$tag,$mfe,$sequence,$fold) = @{$foldInfo};
                my $location = "$chrom:$start..$stop:$strand";
                print ">$tag\t$location\n$sequence\n$fold ($mfe)\n";
            }
        }
    }
}

sub printUniqueHairpinList {
    my $hairpinList = shift;
    foreach my $chrom (keys %{$hairpinList}) {
        foreach my $strand (keys %{$hairpinList->{$chrom}}) {
            # sort by total reads.
	    my @sortedHairpinList = sort {$a->[5] <=> $b->[5]} @{$hairpinList->{$chrom}{$strand}};
            my @uniqueList;
            foreach my $foldInfo (@sortedHairpinList) {
		my($tag,$start,$stop,$hStart,$hStop,$total,$mfe,
		   $aapd,$tapd,$urf,$ahc,$afh,$valid,$shifted) = @{$foldInfo};               
                my $absStart = $start + $hStart;
                my $absStop = $start + $hStop;
                if($strand eq "-") {
                    $absStart = $stop - $hStop;
                    $absStop = $stop - $hStart;
                }
                unless(insideList($absStart,$absStop,\@uniqueList)) {
		    insertEntry(\@uniqueList,[$absStart,$absStop,$foldInfo]);
                }
            }
            foreach my $uniqueFold (@uniqueList) {
                my($absStart,$absStop,$foldInfo) = @{$uniqueFold};
		my($tag,$start,$stop,$hStart,$hStop,$total,$mfe,
		   $aapd,$tapd,$urf,$ahc,$afh,$valid,$shifted) = @{$foldInfo};  
		print "$tag\t$chrom\t$start\t$stop\t$strand\t$hStart\t$hStop\t$total\t$mfe\t$aapd\t$tapd\t$urf\t$ahc\t$afh\t$valid\t$shifted\n";
            }
        }
    }
}

sub extractUniqueHairpinList {
    my $hairpinList = shift;
    my %uniqueHairpinList;
    foreach my $chrom (keys %{$hairpinList}) {
	# sort by total reads, greatest to lowest.
	my @sortedHairpinList = sort {$b->[6] <=> $a->[6]} @{$hairpinList->{$chrom}};
	my @uniqueList;
	foreach my $foldInfo (@sortedHairpinList) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,
	       $aapd,$tapd,$urf,$ahc,$afh,$valid,$shifted) = @{$foldInfo};               
	    my $absStart = $start + $hStart;
	    my $absStop = $start + $hStop;
	    if($strand eq "-") {
		$absStart = $stop - $hStop;
		$absStop = $stop - $hStart;
	    }
	    unless(significantlyInsideList($absStart,$absStop,\@uniqueList)) {
		insertEntry(\@uniqueList,[$absStart,$absStop,$foldInfo]);	
	    }
	}
	foreach my $uniqueFold (@uniqueList) {
	    my($absStart,$absStop,$foldInfo) = @{$uniqueFold};
	    push(@{$uniqueHairpinList{$chrom}},$foldInfo);
	}
    }
    return \%uniqueHairpinList
}

sub multipleHairpins {
    my($fold,$parameters) = @_;
    my $minLength = $parameters->{minLength};
    my $INRUN=0;
    my $lastLeft = 0;
    my @centers;
    for(my $i=0;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq "(") {
	    $lastLeft = $i;
            $INRUN=1;
        } elsif(substr($fold,$i,1) eq ")") {
            if($INRUN) {
                my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$lastLeft,$i);
		my $leftLength = $lastLeft - $leftEnd + 1;
		my $rightLength = $rightEnd - $i + 1;
		#print "$leftLength $rightLength vs. $minLength\n";
		my $length = max($leftLength,$rightLength);
		if($length >= $minLength) {
		    push(@centers,[$lastLeft,$i]);
		}
            }
            $INRUN=0;
	}
    }
    if(scalar(@centers) > 1) {
	return 1;
    }
    return 0;
}

sub getValidHairpinCenters {
    my($fold,$parameters) = @_;
    my $minLength = $parameters->{minLength};
    my $INRUN=0;
    my $lastLeft = 0;
    my @centers;
    for(my $i=0;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq "(") {
            $lastLeft = $i;
            $INRUN=1;
        } elsif(substr($fold,$i,1) eq ")") {
            if($INRUN) {
                my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$lastLeft,$i);
                my $leftLength = $lastLeft - $leftEnd + 1;
                my $rightLength = $rightEnd - $i + 1;
                my $length = max($leftLength,$rightLength);
                if($length > $minLength) {
                    push(@centers,[$lastLeft,$i]);
                }
            }
            $INRUN=0;
        }
    }
    return @centers;
}

sub getHairpinCenters {
    my($fold,$parameters) = @_;
    my $minLength = $parameters->{minLength};
    my $INRUN=0;
    my $lastLeft = 0;
    my @centers;
    for(my $i=0;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq "(") {
            $lastLeft = $i;
            $INRUN=1;
        } elsif(substr($fold,$i,1) eq ")") {
            if($INRUN) {
		push(@centers,[$lastLeft,$i]);
	    }
            $INRUN=0;
        }
    }
    return @centers;
}

sub getMergedHairpinCenters {
    my($fold,$parameters) = @_;
    my $minLength = $parameters->{minLength};
    my $bpDensityLimit = $parameters->{bpDensityLimit};
    my $shortLength = int($minLength*$bpDensityLimit);
    my @centers = getHairpinCenters($fold,$parameters);
    my($mergedCenters,$MERGED) = mergeCenters(\@centers,$fold,$parameters);
    while($MERGED) {
	($mergedCenters,$MERGED) = mergeCenters($mergedCenters,$fold,$parameters);
    }
    return @{$mergedCenters};
}

sub mergeCenters {
    my($centers,$fold,$parameters) = @_;
    my $minLength = $parameters->{minLength};
    my $bpDensityLimit = $parameters->{bpDensityLimit};
    my $shortLength = int($minLength*$bpDensityLimit);
    my @mergedCenters;
    my $MERGED = 0;
    for(my $i=0;$i<@{$centers};$i++) {
	my($leftCenter,$rightCenter) = @{$centers->[$i]};
	#print "checking $leftCenter,$rightCenter\n";
	my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$leftCenter,$rightCenter);
	my $leftLength = $leftCenter - $leftEnd + 1;
	my $rightLength = $rightEnd - $rightCenter + 1;
	if(($leftLength >= $minLength)&&($rightLength >= $shortLength)) {
	    # sufficiently large of a loop, keep as is.
	    push(@mergedCenters,$centers->[$i]);
	} elsif(($leftLength >= $minLength)&&($rightLength < $shortLength)) {
	    # potentially merge the loop.
	    if($i<@{$centers}-1) {
		# if there is a neighboring center
		my($leftCenter2,$rightCenter2) = @{$centers->[$i+1]};
		my($leftEnd2,$rightEnd2) = extendHairpinToFullEnds($fold,$leftCenter2,$rightCenter2);
		my $leftLength2 = $leftCenter2 - $leftEnd2 + 1;
		my $rightLength2 = $rightEnd2 - $rightCenter2 + 1;
		if($leftLength2 < $shortLength) {
		    #print "merging $leftCenter,$rightCenter $leftCenter2,$rightCenter2\n";
		    my($tightLeftEnd,$tightRightEnd) = extendHairpinToEnds($fold,$leftCenter,$rightCenter);
		    my($tightLeftEnd2,$tightRightEnd2) = extendHairpinToEnds($fold,$leftCenter2,$rightCenter2);
		    my $newLeftCenter = getNextLeftPar($fold,$tightLeftEnd);
		    my $newRightCenter = getNextRightPar($fold,$tightRightEnd2);
		    if(($newLeftCenter != -1)&&($newRightCenter != -1)) {
			push(@mergedCenters,[$newLeftCenter,$newRightCenter]);
			$MERGED = 1;
		    }
		}
	    } else {
		push(@mergedCenters,$centers->[$i]);
	    }
	} elsif(($leftLength >= $shortLength)&&($rightLength >= $minLength)) {
	    # sufficiently large of a loop, keep as is.
	    push(@mergedCenters,$centers->[$i]);
	} elsif(($leftLength < $shortLength)&&($rightLength >= $minLength)) {
	    # potentially merge the loop.
	    if($i > 0) {
		# if there is a neighboring center
		my($leftCenter2,$rightCenter2) = @{$centers->[$i-1]};
		my($leftEnd2,$rightEnd2) = extendHairpinToFullEnds($fold,$leftCenter2,$rightCenter2);
		my $leftLength2 = $leftCenter2 - $leftEnd2 + 1;
		my $rightLength2 = $rightEnd2 - $rightCenter2 + 1;
		if($rightLength2 < $shortLength) {
		    #print "merging $leftCenter,$rightCenter $leftCenter2,$rightCenter2\n";
		    my($tightLeftEnd,$tightRightEnd) = extendHairpinToEnds($fold,$leftCenter,$rightCenter);
		    my($tightLeftEnd2,$tightRightEnd2) = extendHairpinToEnds($fold,$leftCenter2,$rightCenter2);
		    my $newLeftCenter = getNextLeftPar($fold,$tightLeftEnd2);
		    my $newRightCenter = getNextRightPar($fold,$tightRightEnd);
		    if(($newLeftCenter != -1)&&($newRightCenter != -1)) {
			push(@mergedCenters,[$newLeftCenter,$newRightCenter]);
			$MERGED = 1;
		    }
		}
	    } else {
		push(@mergedCenters,$centers->[$i]);
	    }
	}
    }
    my @filteredCenters;
    for(my $i=0;$i<@mergedCenters;$i++) {
	my($l1,$r1) = @{$mergedCenters[$i]};
	my $N = @filteredCenters;
	my $REMOVE=0;
	for(my $j=0;$j<$N;$j++) {
	    my($l2,$r2) = @{$filteredCenters[$j]};
	    if(($l1==$l2)&&($r1==$r2)) {
		$REMOVE=1;
	    }
	}
	unless($REMOVE) {
	    push(@filteredCenters,[$l1,$r1]);
	}
    }
    return(\@filteredCenters,$MERGED);
}

sub getMirHairpinCenters {
    my($fold,$parameters) = @_;
    my $minLength = $parameters->{minLength};
    my $INRUN=0;
    my $lastLeft = 0;
    my @centers;
    for(my $i=0;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq "(") {
            $lastLeft = $i;
            $INRUN=1;
        } elsif(substr($fold,$i,1) eq ")") {
            if($INRUN) {
		push(@centers,[$lastLeft,$i]);
            }
            $INRUN=0;
        }
    }
    if(@centers > 1) {
	my @newCenters;
	for(my $i=0;$i<@centers-1;$i++) {
	    my($left1,$right1) = @{$centers[$i]};
	    my($left2,$right2) = @{$centers[$i+1]};
	    # check the separation
	    if($left2 - $right1 < 40) {
		my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$left1,$right2);
                my $leftLength = $left1 - $leftEnd + 1;
                my $rightLength = $rightEnd - $right2 + 1;
                my $length = max($leftLength,$rightLength);
                if($length > $minLength) {
		    push(@newCenters,[$left1,$right2]);
		    $i++;
		}
	    } else {
		my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$left1,$right1);
                my $leftLength = $left1 - $leftEnd + 1;
                my $rightLength = $rightEnd - $right1 + 1;
                my $length = max($leftLength,$rightLength);
                if($length > $minLength) {
		    push(@newCenters,[$left1,$right1]);
		}
	    }
	}
	return @newCenters;
    } 
    return @centers;
}

sub getFullHairpinWindow {
    # this method returns the window of paired bases, extends beyond minor hairpins    
    my($fold,$leftCenter,$rightCenter) = @_;;
    return extendHairpinToFullEnds($fold,$leftCenter,$rightCenter);
}

sub getHairpinWindow {
    # this method returns the window of paired bases, but stops at minor hairpins.
    my($fold,$leftCenter,$rightCenter) = @_;
    return extendHairpinToEnds($fold,$leftCenter,$rightCenter);
}

sub getMaxLeftCenter {
    my $fold = shift;
    my $INRUN=0;
    my $maxLength = 0;
    my $maxRunPos=0;
    my $lastLeft = 0;
    for(my $i=0;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq "(") {
	    $lastLeft = $i;
            $INRUN=1;
        } elsif(substr($fold,$i,1) eq ")") {
            if($INRUN) {
                my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$lastLeft,$i);
                my $length = $rightEnd - $leftEnd + 1;
                if($length > $maxLength) {
                    $maxLength = $length;
                    $maxRunPos = $lastLeft;
                }
            }
            $INRUN=0;
	}
    }
    return $maxRunPos;
}

sub getMaxRightCenter {
    my $fold = shift;
    my $leftCenter = getMaxLeftCenter($fold);
    for(my $i=$leftCenter;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq ")") {
            return $i;
        }
    }
    return length($fold)-1;
}

sub getFullBasePairMap {
    my($leftCenter,$rightCenter,$sequence,$parameters) = @_;
    my $RNAfold = $parameters->{RNAfold};
    my $fastaFile = "temp.fasta";                                 
    my $foldFile = "temp.mfe";
    my $dotFile = "temp_dp.ps";
    # print the sequence for folding.
    open(TF,">$fastaFile");
    print TF ">temp\n$sequence\n@\n";
    close(TF);
    # fold the UTR context region.
    system("cat $fastaFile | $RNAfold -noPS -p > $foldFile");
    # read the energies of the folds.
    my @basePairMap;
    open(DF,"$dotFile");
    while(<DF>) {
	chomp;
	next unless /(\d+) (\d+) (\d+\.\d+) (.box)$/;
        my ($i, $j, $p, $id) = ($1,$2,$3,$4);
	if($id eq "lbox") {
	    $i = $i-1;
	    $j = $j-1;
	    push(@basePairMap,[$i,$j]);
	}	
    }
    close(DF);
    # reverse the order read in...
    @basePairMap = reverse @basePairMap;
    return \@basePairMap;
}

sub getBasePairMap {
    my($leftCenter,$rightCenter,$sequence,$parameters) = @_;
    my $RNAfold = $parameters->{RNAfold};
    my $fastaFile = "temp.fasta";                                 
    my $foldFile = "temp.mfe";
    my $dotFile = "temp_dp.ps";
    # print the sequence for folding.
    open(TF,">$fastaFile");
    print TF ">temp\n$sequence\n@\n";
    close(TF);
    # fold the UTR context region.
    system("cat $fastaFile | $RNAfold -noPS -p > $foldFile");
    # read the energies of the folds.
    my @basePairMap;
    open(DF,"$dotFile");
    while(<DF>) {
	chomp;
	next unless /(\d+) (\d+) (\d+\.\d+) (.box)$/;
        my ($i, $j, $p, $id) = ($1,$2,$3,$4);
	if($id eq "lbox") {
	    $i = $i-1;
	    $j = $j-1;
	    if(($i<=$leftCenter)&&($rightCenter<=$j)) {
		push(@basePairMap,[$i,$j]);
	    }	
	}
    }
    close(DF);
    # reverse the order read in...
    @basePairMap = reverse @basePairMap;
    return \@basePairMap;
}

sub getBasePairs {
    my($leftCenter,$rightCenter,$fold) = @_;
    my @basePairs;
    my $leftCurrent = $leftCenter+1;
    my $rightCurrent = $rightCenter-1;
    my $STAY=1;
    my $COUNT=0;
    while($STAY) {
        my $left = getNextLeftPar($fold,$leftCurrent);
        my $right = getNextRightPar($fold,$rightCurrent);
        if(($left == -1)||($right == -1)) {
            $STAY=0;
        } else {
            push(@basePairs,[$left,$right]);
            $COUNT++;
            $leftCurrent = $left;
            $rightCurrent = $right;
        }
    }
    if(@basePairs) {
        return \@basePairs;
    } else {
        die "crapped out at $leftCenter,$rightCenter:\n$fold\n";
    }
}

sub getMaxHairpinLength {
    my($fold,$leftCenter,$rightCenter) = @_;
    my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$leftCenter,$rightCenter);
    my $leftLength = $leftCenter - $leftEnd + 1;
    my $rightLength = $rightEnd - $rightCenter + 1;
    return max($leftLength,$rightLength);
}

sub extendHairpinToEnds {
    my($fold,$leftCenter,$rightCenter) = @_;
    my @basePairs;
    # start left current one up, so the true left middle will be found first.
    my $leftCurrent = $leftCenter+1;
    # start the search for the right par on the next base.
    my $rightCurrent = $rightCenter-1;
    my $STAY=1;
    my $COUNT=0;
    while($STAY) {
        my $left = getNextLeftPar($fold,$leftCurrent);
        my $right = getNextRightPar($fold,$rightCurrent);
        if(($left == -1)||($right == -1)) {
            $STAY=0;
	} else {
            push(@basePairs,[$left,$right]);
	    $COUNT++;
            $leftCurrent = $left;
            $rightCurrent = $right;
        }
    }
    unless(@basePairs) {
	my $L = "L" x ($leftCenter+1);
	my $R = "R" x ($rightCenter+1);
	die "extendHairpinToEnds: FAIL:\n$fold\t$L\n$R\n";
    }
    my $endPairs = pop(@basePairs);
    return @{$endPairs};
}

sub extendHairpinToFullEnds {
    my($fold,$leftCenter,$rightCenter) = @_;
    my @basePairs;
    # start left current one up, so the true left middle will be found first.  
    my $leftCurrent = $leftCenter+1;
    # start the search for the right par on the next base.
    my $rightCurrent = $rightCenter-1;
    my $STAY=1;
    my $COUNT=0;
    while($STAY) {
        my $left = getNextLeftPar($fold,$leftCurrent);
        my $right = getNextRightPar($fold,$rightCurrent);
        if(($left == -1)||($right == -1)) {
            $STAY=0;
        } else {
            push(@basePairs,[$left,$right]);
	    $COUNT++;
            $leftCurrent = $left;
            $rightCurrent = $right;
        }
    }
    unless(@basePairs) {
	my $L = "L" x ($leftCenter+1);
	my $R = "R" x ($rightCenter+1);
        die "extendHairpinFullToEnds FAIL: ($leftCenter,$rightCenter)\n$fold\n$L\n$R\n\n";	
    }
    my $endPairs = pop(@basePairs);
    my($leftEnd,$rightEnd) = @{$endPairs};
    while(($leftEnd > 0)&&(substr($fold,$leftEnd-1,1) ne ")")) {
	$leftEnd--;
    }
    while(($rightEnd < length($fold)-1)&&(substr($fold,$rightEnd+1,1) ne "(")) {
	$rightEnd++;
    }
    while(substr($fold,$leftEnd,1) eq ".") {
	$leftEnd++;
    }
    while(substr($fold,$rightEnd,1) eq ".") {
	$rightEnd--;
    }
    return ($leftEnd,$rightEnd);
}

sub getNextLeftPar {
    my($fold,$start) = @_;
    for(my $i=$start-1;$i>=0;$i--) {
        if(substr($fold,$i,1) eq "(") {
            return $i;
        } elsif(substr($fold,$i,1) eq ")") {
            return -1;
        }
    }
    return -1;
}

sub getNextRightPar {
    my($fold,$start)=@_;
    for(my $i=$start+1;$i<length($fold);$i++) {
	if(substr($fold,$i,1) eq ")") {
            return $i;
        } elsif(substr($fold,$i,1) eq "(") {
	    return -1
	}
    }
    return -1;
}

###############################
# READ PROCESSING SUBROUTINES #
###############################

sub printAllDistinctReads {
    my($dataSet,$sampleList,$hitCount,$parameters) = @_;
    open(DR,">distinctReads.txt");
    print DR "#chrom\tstrand\tstart\tstop\ttotal\ttotalHitCount\treadTotal";
    foreach my $sample (@{$sampleList}) {
	print DR "\t$sample";
    }
    print DR "\n";
    foreach my $chrom (keys %{$dataSet}) {
	foreach my $strand (keys %{$dataSet->{$chrom}}) {
	    my %distinctReadCount;
	    my %libraryCounts;
	    my %totalHitCount;
	    my %distinctReadTotal;
	    foreach my $read (@{$dataSet->{$chrom}{$strand}}) {
		my($rStart,$rStop,$sample,$rID) = @{$read};
		$distinctReadCount{$rStart."x".$rStop} += 1/($hitCount->{$sample.$rID});
		$totalHitCount{$rStart."x".$rStop} += $hitCount->{$sample.$rID};
		$distinctReadTotal{$rStart."x".$rStop}++;
		$libraryCounts{$rStart."x".$rStop}{$sample} += 1/($hitCount->{$sample.$rID});
	    }
	    my @distinctReads;
	    foreach my $key (keys %distinctReadCount) {
		my($drStart,$drStop) = split(/x/,$key);
		foreach my $sample (@{$sampleList}) {
		    unless($libraryCounts{$key}{$sample}) {
			$libraryCounts{$key}{$sample} = 0;
		    }
		}		
		push(@distinctReads,[$drStart,$drStop,$distinctReadCount{$key},$libraryCounts{$key},$totalHitCount{$key},$distinctReadTotal{$key}]);
	    }
	    @distinctReads = sort {$a->[0] <=> $b->[0]} @distinctReads;
	    foreach my $distinctRead (@distinctReads) {
		my($drStart,$drStop,$total,$libraryCounts,$totalHitCount,$readTotal) = @{$distinctRead};
		print DR "$chrom\t$strand\t$drStart\t$drStop\t";
		printf(DR "%.3f",$total);
		print DR "\t$totalHitCount\t$readTotal";
		foreach my $sample (@{$sampleList}) {
		    printf(DR "\t%.3f", $libraryCounts->{$sample});
		}
		print DR "\n";
	    }
	}
    }
    close(DR);
}

sub processOverlaps {
    my($dataSet,$repeatRegions,$hitCount,$parameters)=@_;
    my $maxLength = $parameters->{maxLength} or die "FAIL: maxLength not loaded.\n";
    my $distancesFile = $parameters->{distancesFile} or die "FAIL: distancesFile not loaded.\n";
    my $readRegionsFile = $parameters->{readRegionsFile} or die "FAIL: readRegionsFile not loaded.\n";
    my $longReadRegionsFile = $parameters->{longReadRegionsFile} or die "FAIL: longReadRegionsFile not loaded.\n";
    my $allReadRegionsFile = $parameters->{allReadRegionsFile} or die "FAIL: allReadRegionsFile not loaded.\n";
    open(DF,">".$distancesFile);
    open(RF,">".$readRegionsFile);
    open(LRF,">".$longReadRegionsFile);
    open(ARF,">".$allReadRegionsFile);
    # for regions shorter than 100bp.
    my $COUNT1 = 1;
    # for regions less than 100bp.
    my $COUNT2 = 1;
    # for repeat regions.
    my $COUNT3++;
    foreach my $chr (keys %{$dataSet}) {
        foreach my $strand (keys %{$dataSet->{$chr}}) {
            # extract regions in a strand-dependent way.
            my $overlapList = extractOverlapList($chr,$strand,
						 $dataSet->{$chr}{$strand},$hitCount,$parameters);
            # print overlap regions.
            my $prevStop;
            # for each region.
            foreach my $region (@{$overlapList}) {
                my($start,$stop) = @{$region};
                # unless it overlaps a repeat region.
		my $location = "$chr:$start..$stop:$strand";
		my $length = $stop - $start + 1;
		unless(insideList($start,$stop,$repeatRegions->{$chr})) {
                    # print the info.
                    if($length <= $maxLength) {
			print ARF "mmr$COUNT1\t$location\t$length\n";
                        print RF "mmr$COUNT1\t$location\t$length\n";
                        $COUNT1++;
                    } else {
			print ARF "lrr$COUNT2\t$location\t$length\trejected: region too long.\n";
                        print LRF "lrr$COUNT2\t$location\t$length\n";
                        $COUNT2++;
                    }
                    if($prevStop) {
                        my $dist = $start - $prevStop + 1;
                        print DF "$dist\n";
                    }
                    $prevStop = $stop;
                } else {		    
		    print ARF "rrr$COUNT3\t$location\t$length\trejected: region within repeat.\n";
		    $COUNT3++;
		}
            }
        }
    }
    close(DF);
    close(RF);
    close(LRF);
    close(ARF);
}

sub extractOverlapList {
    # this is run for a particular chrom and strand
    # chromDataSet is for a particular chrom and strand of the readDataSet
    my($chr,$strand,$chromDataSet,$hitCount,$parameters)=@_;
    my $maxCount = $parameters->{maxCount};
    my @overlapList;
    my $maxJ = 0; # to store the highest position covered by a read on the chrom/strand
    my @overlaps;
    foreach my $site (@{$chromDataSet}) {
        my($rStart,$rStop,$sample,$id)=@{$site};
        if($hitCount->{$sample.$id} <= $maxCount) {
            for(my $j=$rStart;$j<=$rStop;$j++) {
                $overlaps[$j]++; # here is the increment for a given read.
                if($maxJ < $j) {
                    $maxJ = $j;
                }
            }
        }
    }
    $maxJ += 100;
    my $prevStart = 1;
    my $prevStop = 1;
    for(my $j=1;$j<=$maxJ;$j++) {
        if($overlaps[$j]) {
            # if still at non-zero value, extend the window.
            $prevStop = $j;
        } else {
            # if at zero value, store the previous window information, and start a new one.
            if($prevStart != $prevStop) {
                push(@overlapList,[$prevStart,$prevStop]);
            }
            $prevStart = $j;
            $prevStop = $j; # for now, assume the overlap is only one base.
        }
    }
    return \@overlapList;
}

sub printReadRegionRNAfold {
    my($readRegionsFile,$genome,$parameters)=@_;
    my $readRegionsFastaFile = $parameters->{readRegionsFastaFile} or die "readRegionsFastaFile not loaded.\n";
    my $readRegions = readReadRegions($readRegionsFile);
    open(RNAF,">$readRegionsFastaFile") or die "could not open $readRegionsFastaFile for writing.\n";
    foreach my $region (@{$readRegions}) {
        my($id,$location,$length) = @{$region};
        my($chrom,$start,$stop,$strand) = parseLocation($location);
        my $bufferLength = getBufferLength($start,$stop,$parameters);
        my $newStart = $start-$bufferLength > 1 ? $start - $bufferLength : 1;
        my $newStop = $stop+$bufferLength <= length($genome->{$chrom}) ? 
	    $stop+$bufferLength : 
	    length($genome->{$chrom});
        my $newLocation = $chrom . ":" . $newStart . ".." . $newStop  . ":" . $strand;
        my $sequence = getSequence($newLocation,$genome);
        print RNAF ">$id\t$newLocation\n$sequence\n";
    }
    print RNAF "@\n";
    close(RNAF);
}

sub printMirRegionRNAfold {
    my($mirRegionsFile,$genome,$parameters)=@_;
    my $mirRegionsFastaFile = $parameters->{mirRegionsFastaFile} or die "mirRegionsFastaFile not defined.\n";
    my $mirRegions = readReadRegions($mirRegionsFile);
    my $buffer = 20;
    open(RNAF,">$mirRegionsFastaFile") or die "could not open $mirRegionsFastaFile for writing.\n";
    foreach my $region (@{$mirRegions}) {
        my($id,$location,$length) = @{$region};
	my($chrom,$start,$stop,$strand) = parseLocation($location);
	my $newLocation = $chrom . ":" . ($start-$buffer) . ".." . ($stop+$buffer)  . ":" . $strand;
        my $sequence = getSequence($newLocation,$genome);
        print RNAF ">$id\t$location\n$sequence\n";
    }
    print RNAF "@\n";
    close(RNAF);
}

sub printMirLociToMirList {
    my($readRegionsFile,$genome,$parameters)=@_;
    my $readRegions = readReadRegions($readRegionsFile);
    my $RNAfold = $parameters->{RNAfold};
    open(RNAF,"$readRegionsFile");
    while(<RNAF>) {
        my($id,$location,$strand,$total) = split(/\t/);;
	my($chrom,$start,$stop) = parseLocation($location);
	my $newLocation = $location . ":" . $strand;
        my $sequence = getSequence($newLocation,$genome);
	my $fastaFile = "fasta.fasta";                                 
	my $foldFile = "fold.mfe";
	open(TF,">$fastaFile");
	print TF ">temp\n$sequence\n@\n";
	close(TF);
	# fold the UTR context region.
	system("cat $fastaFile | $RNAfold -noPS > $foldFile");
	# read the energies of the folds.
	my($fold,$mfe);
	open(TMFE,"$foldFile");
	while(<TMFE>) {
	    chomp;
	    if(/[ACGU]/) {
		my $seq = $_;
		$seq =~ tr/acgtuU/ACGTTT/;
		if($seq ne $sequence) {
		    die "error reading RNAfold output.\n$seq != $sequence\n";
		}
	    } elsif(/(.*)\s+\((\s*-*\d+\.\d+)\)/) {
		$fold = $1;
		$mfe = $2;           
	    }
	}
	close(TMFE);
	print "$id\t$newLocation\t";
    }
    print RNAF "@\n";
    close(RNAF);
}

sub getBufferLength {
    my($start,$stop,$parameters) = @_;
    my $totalLength = $parameters->{totalLength};
    if($stop - $start + 1 > 50 ) {
	if($stop - $start + 1 < $totalLength) {
	    return int(($totalLength - ($stop - $start + 1))/2);
	} else {
	    return 0;
	}
    } else {
        return int($totalLength/2);
    }
}

sub processReadRegions {
    my($folds,$dataSet,$sampleList,$hitCount,$tRNAs,$parameters)=@_;
    my $minLength = $parameters->{minLength} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile} or die "filteredCandidatesFile: not loaded.\n";
    my $hairpinsFile = $parameters->{hairpinsFile} or die "hairpinsFile: not loaded.\n";
    my $distinctReadsFile = $parameters->{distinctReadsFile} or die "distinctReadsFile: not loaded.\n";
    my $productFile = $parameters->{productFile} or die "productFile: not loaded.\n";
    open(FRC,">".$filteredCandidatesFile);
    close(FRC);
    open(HL,">".$hairpinsFile);
    open(HDR,">".$distinctReadsFile);
    open(HPL,">".$productFile);
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalReads\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tafh\tvalid\tsameShifted\tbothShift\n";
    print HDR "#tag\tstart\tstop\ttotal reads\toffsetStart\trelative Start";
    foreach my $sample (@{$sampleList}) {
	print HDR "\t$sample";
    }
    print HDR "\n";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tsequence";
    foreach my $sample (@{$sampleList}) {
        print HPL "\t$sample";
    }
    print HPL "\n";
    $tRNAs = processTrnaOverlaps($folds,$tRNAs);
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand (keys %{$folds->{$chrom}}) {	    
	    foreach my $foldInfo (@{$folds->{$chrom}{$strand}}) {
		my($start,$stop,$tag,$mfe,$sequence,$fold) = @{$foldInfo};
		my $location = "$chrom:$start..$stop:$strand";
		#print "checking $tag\n";
		my @centers = getMergedHairpinCenters($fold,$parameters);
		my @labels = ("a","b","c","d","e","f","g","h","i");
		my $COUNT=0;
		foreach my $center (@centers) {
		    my($leftCenter,$rightCenter) = @{$center};
		    #print "processing $leftCenter,$rightCenter\n";
		    my $newId = $tag . $labels[$COUNT];
		    $COUNT++;
		    if(goodCandidate($tag,$location,$tRNAs,$parameters)) {
			#my($chrom,$start,$stop,$strand) = parseLocation($location);
			#print "$tag.$labels[$COUNT]\n";
			my $basePairs = getBasePairs($leftCenter,$rightCenter,
						     $fold,$parameters);
			my $hairpinLength = getMaxHairpinLength($fold,$leftCenter,$rightCenter);
			my($hStart,$hStop) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
			#print "hp: $maxHairpinLength $hStart $hStop $minLength\n";
			if($hairpinLength >= $minLength) {
			    my $distinctReads=getDistinctReadInfo($start,$stop,
								  $dataSet->{$chrom}{$strand},
								  $sampleList,$hitCount);
			    my $productInfo=extractProducts($leftCenter,$rightCenter,
							    $fold,$basePairs,$location,
							    $distinctReads,$parameters);
			    if(goodProducts($productInfo,$parameters)) {							    
				my $revStrand = revStrand($strand);
				my $revDistinctReads=getDistinctReadInfo($start,$stop,
									 $dataSet->{$chrom}{$revStrand},
									 $sampleList,$hitCount);
				my $revProductInfo=extractProducts($leftCenter,$rightCenter,
								   $fold,$basePairs,$location,
								   $revDistinctReads,$parameters);
				#doubleStranded($start,$stop,$strand,$dataSet->{$chrom},$hitCount,$parameters);
				my($PASS,$REASON) = plausibleReads($productInfo,$revProductInfo,$parameters);
				if($PASS) {
				    my($tpd,$totalRP)=getReverseProductDisplacement($productInfo,
										    $revProductInfo,
										    $parameters);
				    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
				    my $urf = computeUniqueReadFraction($tag,$location,
									$dataSet,$hitCount);
				    my $ahc = computeMaxProdHitCount($productInfo,$location,
								     $dataSet,$hitCount,$parameters);
				    my $afh = computeMaxProdFivePrimeHet($productInfo,$parameters);
				    my $pbp = computeProductBasePairing($leftCenter,$rightCenter,
									$productInfo,$basePairs,$parameters);
				    my $sameShift = computeMaxSameShift($leftCenter,$rightCenter,
									$fold,$basePairs,
									$location,$distinctReads,$parameters);
				    my $bothShift = computeMaxBothShift($leftCenter,$rightCenter,
									$fold,$basePairs,
									$location,$distinctReads,$parameters);
				    # total is the read count for the whole hairpin
				    my $total = 0;
				    foreach my $product (@{$productInfo}) {
					my($side,$newType,$prodList,$prodCount,$maxProdCount,
					   $relStart,$relStop,$offset,$gStart) = @{$product};
					my $length = $relStop - $relStart + 1;
					my $productSequence = substr($sequence,$relStart,$length);
					$total += $prodCount;
					my $totalLibraryCounts = getTotalLibraryCounts($prodList);
					print HPL "$newId\t$side\t$newType\t$prodCount\t";
					print HPL "$maxProdCount\t$relStart\t$relStop\t";
					print HPL "$productSequence";				
					foreach my $sample (@{$sampleList}) {
					    printf(HPL "\t%.3f",$totalLibraryCounts->{$sample});
					}
					print HPL "\n";
					foreach my $read (@{$prodList}) {
					    my($relStart,$relStop,$offset,$gStart,
					       $count,$libraryCounts) = @{$read};
					    my $gStop = $gStart + ($relStop - $relStart);
					    print HDR "$newId\t$gStart\t$gStop\t$count\t";
					    print HDR "$offset\t$relStart";
					    foreach my $sample (@{$sampleList}) {
						printf(HDR "\t%.3f",$libraryCounts->{$sample});
					    }
					    print HDR "\n";
					}
				    }
				    print HL "$newId\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
				    printf(HL "%.3f\t%.3f\t",$total,$mfe);
				    print HL "$sequence\t$fold\t";
				    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
					   $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift);
				} else {
				    open(FRC,">>$parameters->{filteredCandidatesFile}");
				    print FRC "$newId\t$location\t$REASON\n";
				    close(FRC);
				}
			    } else {
				open(FRC,">>$parameters->{filteredCandidatesFile}");
				print FRC "$newId\t$location\trejected: no good products.\n";
				close(FRC);
			    }
			} else {
			    open(FRC,">>$parameters->{filteredCandidatesFile}");
			    print FRC "$newId\t$location\trejected: longest hairpin too short. ";
			    print FRC "max length=$hairpinLength\n";
			    close(FRC);
			}
		    }
		}
	    }
	}
    }
}

sub processReadRegions2 {
    my($folds,$allDistinctReads,$sampleList,$tRNAs,$parameters)=@_;   
    my $minLength = $parameters->{minLength};
    open(FRC,">$parameters->{filteredCandidatesFile}");
    close(FRC);
    open(HL,">hairpinList.txt");
    open(HDR,">hairpinDistinctReads.txt");
    open(HPL,">hairpinProductList.txt");
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalReads\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tafh\tvalid\tsameShift\tbothShift\n";
    print HDR "#tag\tstart\tstop\ttotal reads\toffsetStart\trelative Start";
    foreach my $sample (@{$sampleList}) {
	print HDR "\t$sample";
    }
    print HDR "\n";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tsequence";
    foreach my $sample (@{$sampleList}) {
        print HPL "\t$sample";
    }
    print HPL "\n";
    $tRNAs = processTrnaOverlaps($folds,$tRNAs);
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand (keys %{$folds->{$chrom}}) {	    
	    foreach my $foldInfo (@{$folds->{$chrom}{$strand}}) {
		my($start,$stop,$tag,$mfe,$sequence,$fold) = @{$foldInfo};
		my $location = "$chrom:$start..$stop:$strand";
		my @centers = getMergedHairpinCenters($fold,$parameters);
		my @labels = ("a","b","c","d","e","f","g","h");
		my $COUNT=0;
		foreach my $center (@centers) {
		    my($leftCenter,$rightCenter) = @{$center};
		    my $newId = $tag . $labels[$COUNT];
		    $COUNT++;
		    if(goodCandidate($tag,$location,$tRNAs,$parameters)) {
			#my($chrom,$start,$stop,$strand) = parseLocation($location);
			#print "$tag.$labels[$COUNT]\n";
			my $basePairs = getBasePairs($leftCenter,$rightCenter,
						     $fold,$parameters);
			my $hairpinLength = getMaxHairpinLength($fold,$leftCenter,$rightCenter);
			my($hStart,$hStop) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
			#print "hp: $maxHairpinLength $hStart $hStop $minLength\n";
			if($hairpinLength >= $minLength) {
			    my $distinctReads=extractDistinctReadInfo($start,$stop,
								      $allDistinctReads->{$chrom}{$strand});			    
			    my $productInfo=extractProducts($leftCenter,$rightCenter,
							    $fold,$basePairs,$location,
							    $distinctReads,$parameters);
			    if(goodProducts($productInfo,$parameters)) {
				my $revStrand = revStrand($strand);
				my $revDistinctReads=extractDistinctReadInfo($start,$stop,
									     $allDistinctReads->{$chrom}{$revStrand});
				my $revProductInfo=extractProducts($leftCenter,$rightCenter,
								   $fold,$basePairs,$location,
								   $revDistinctReads,$parameters);
				#doubleStranded($start,$stop,$strand,$dataSet->{$chrom},$hitCount,$parameters);
				my($PASS,$REASON) = plausibleReads($productInfo,$revProductInfo,$parameters);
				if($PASS) {
				    my($tpd,$totalRP)=getReverseProductDisplacement($productInfo,
										    $revProductInfo,
										    $parameters);
				    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
				    my $urf = computeUniqueReadFraction2($distinctReads);
				    my $ahc = computeMaxProdHitCount2($productInfo,$distinctReads,$parameters);
				    my $afh = computeMaxProdFivePrimeHet($productInfo,$parameters);
				    my $pbp = computeProductBasePairing($leftCenter,$rightCenter,
									$productInfo,$basePairs,$parameters);
				    my $sameShift = computeMaxSameShift($leftCenter,$rightCenter,
									$fold,$basePairs,
									$location,$distinctReads,$parameters);
				    my $bothShift = computeMaxBothShift($leftCenter,$rightCenter,
									$fold,$basePairs,
									$location,$distinctReads,$parameters);
				    # total is the read count for the whole hairpin
				    my $total = 0;
				    foreach my $product (@{$productInfo}) {
					my($side,$newType,$prodList,$prodCount,$maxProdCount,
					   $relStart,$relStop,$offset,$gStart) = @{$product};
					my $length = $relStop - $relStart + 1;
					my $productSequence = substr($sequence,$relStart,$length);
					$total += $prodCount;
					my $totalLibraryCounts = getTotalLibraryCounts($prodList);
					print HPL "$newId\t$side\t$newType\t$prodCount\t";
					print HPL "$maxProdCount\t$relStart\t$relStop\t";
					print HPL "$productSequence";
					foreach my $sample (@{$sampleList}) {
					    printf(HPL "\t%.3f",$totalLibraryCounts->{$sample});
					}
					print HPL "\n";
					foreach my $read (@{$prodList}) {
					    my($relStart,$relStop,$offset,$gStart,
					       $count,$libraryCounts) = @{$read};
					    my $gStop = $gStart + ($relStop - $relStart);
					    print HDR "$newId\t$gStart\t$gStop\t$count\t";
					    print HDR "$offset\t$relStart";
					    foreach my $sample (@{$sampleList}) {
						printf(HDR "\t%.3f",$libraryCounts->{$sample});
					    }
					    print HDR "\n";
					}
				    }
				    print HL "$newId\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter";
				    printf(HL "%.3f\t%.3f\t",$total,$mfe);
				    print HL "$sequence\t$fold\n";
				    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t",
					   $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift);
				} else {
				    open(FRC,">>$parameters->{filteredCandidatesFile}");
				    print FRC "$newId\t$location\t$REASON\n";
				    close(FRC);
				}
			    } else {
				open(FRC,">>$parameters->{filteredCandidatesFile}");
				print FRC "$newId\t$location\trejected: no good products.\n";
				close(FRC);
			    }
			} else {
			    open(FRC,">>$parameters->{filteredCandidatesFile}");
			    print FRC "$newId\t$location\trejected: longest hairpin too short. ";
			    print FRC "max length=$hairpinLength\n";
			    close(FRC);
			}
		    }
		}
	    }
	}
    }
}

sub processMirRegions {
    my($folds,$dataSet,$sampleList,$hitCount,$parameters)=@_;   
    my $minLength = $parameters->{minLength};
    open(HL,">miRList.txt");
    open(HDR,">miRDistinctReads.txt");
    open(HPL,">miRProductList.txt");
    print HL "#tag\tchrom\tstart\tstop\tstrand\thStart\thStop\t";
    print HL "totalReads\tmfe\tsequence\tfold\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotalReads\toffsetStart\trelative Start";
    foreach my $sample (@{$sampleList}) {
	print HDR "\t$sample";
    }
    print HDR "\n";
    print HPL "#tag\tside\ttype\ttotalReads\ttotalMostAbundant\tstart\tstop\tsequence";
    foreach my $sample (@{$sampleList}) {
        print HPL "\t$sample";
    }
    print HPL "\n";
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand (keys %{$folds->{$chrom}}) {	    
	    foreach my $foldInfo (@{$folds->{$chrom}{$strand}}) {
		my($tag,$start,$stop,$leftCenter,$rightCenter,$mfe,$sequence,$fold) = @{$foldInfo};
		#print "checking $tag ";
		my $location = "$chrom:$start..$stop:$strand";		
		my $basePairs = getBasePairs($leftCenter,$rightCenter,
					     $fold,$parameters);
		my $hairpinLength = getMaxHairpinLength($fold,$leftCenter,$rightCenter);
		my($hStart,$hStop) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
		my $distinctReads=getDistinctReadInfo($start,$stop,
						      $dataSet->{$chrom}{$strand},
						      $sampleList,$hitCount);		
		my $productInfo=extractProducts($leftCenter,$rightCenter,
						$fold,$basePairs,$location,
						$distinctReads,$parameters);
		my $revStrand = revStrand($strand);
		my $revDistinctReads=getDistinctReadInfo($start,$stop,
							 $dataSet->{$chrom}{$revStrand},
							 $sampleList,$hitCount);
		my $revProductInfo=extractProducts($leftCenter,$rightCenter,
						   $fold,$basePairs,$location,
						   $revDistinctReads,$parameters);
		my $total = 0;
		foreach my $product (@{$productInfo}) {
		    my($side,$newType,$prodList,$prodCount,$maxProdCount,
		       $relStart,$relStop,$offset,$gStart) = @{$product};
		    my $length = $relStop - $relStart + 1;
		    my $productSequence = substr($sequence,$relStart,$length);
		    unless($newType eq "out") {
			$total += $prodCount;
		    }
		    my $totalLibraryCounts = getTotalLibraryCounts($prodList);
		    print HPL "$tag\t$side\t$newType\t$prodCount\t";
		    print HPL "$maxProdCount\t$relStart\t$relStop\t";
		    print HPL "$productSequence";				
		    foreach my $sample (@{$sampleList}) {
			printf(HPL "\t%.3f",$totalLibraryCounts->{$sample});
		    }
		    print HPL "\n";
		    foreach my $read (@{$prodList}) {
			my($relStart,$relStop,$offset,$gStart,
			   $count,$libraryCounts) = @{$read};
			my $gStop = $gStart + ($relStop - $relStart);
			print HDR "$tag\t$gStart\t$gStop\t$strand\t$count\t";
			print HDR "$offset\t$relStart";
			foreach my $sample (@{$sampleList}) {
			    printf(HDR "\t%.3f",$libraryCounts->{$sample});
			}
			print HDR "\n";
		    }
		}
		foreach my $product (@{$revProductInfo}) {
		    my($side,$newType,$prodList,$prodCount,$maxProdCount,
		       $relStart,$relStop,$offset,$gStart) = @{$product};
		    $newType .= "roR";
		    my $length = $relStop - $relStart + 1;
		    my $productSequence = substr($sequence,$relStart,$length);
		    my $totalLibraryCounts = getTotalLibraryCounts($prodList);
		    print HPL "$tag\t$side\t$newType\t$prodCount\t";
		    print HPL "$maxProdCount\t$relStart\t$relStop\t";
		    print HPL "$productSequence";				
		    foreach my $sample (@{$sampleList}) {
			printf(HPL "\t%.3f",$totalLibraryCounts->{$sample});
		    }
		    print HPL "\n";
		    foreach my $read (@{$prodList}) {
			my($relStart,$relStop,$offset,$gStart,
			   $count,$libraryCounts) = @{$read};
			my $gStop = $gStart + ($relStop - $relStart);
			print HDR "$tag\t$gStart\t$gStop\t$revStrand\t$count\t";
			print HDR "$offset\t$relStart";
			foreach my $sample (@{$sampleList}) {
			    printf(HDR "\t%.3f",$libraryCounts->{$sample});
			}
			print HDR "\n";
		    }
		}
		print HL "$tag\t$chrom\t$start\t$stop\t$strand\t$hStart\t$hStop\t";
		printf(HL "%.3f\t%.3f\t",$total,$mfe);
		print HL "$sequence\t$fold\n";
	    }
	}
    }		    
}

sub processTrnaOverlaps {
    my($folds,$tRNAs) = @_;
    my %newTrnas;
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand1 (keys %{$folds->{$chrom}}) {
	    foreach my $foldInfo1 (@{$folds->{$chrom}{$strand1}}) {
		my($start1,$stop1,$tag1) = @{$foldInfo1};
		if($tRNAs->{$tag1}) {
		    $newTrnas{$tag1} = 1;
		    foreach my $strand2 (keys %{$folds->{$chrom}}) {	    
			foreach my $foldInfo2 (@{$folds->{$chrom}{$strand2}}) {
			    my($start2,$stop2,$tag2) = @{$foldInfo2};
			    if(getOverlap($start1,$stop1,$start2,$stop2)) {
				$newTrnas{$tag2} = 1;
			    }
			}
		    }
		}
	    }
	}
    }
    return \%newTrnas;
}

sub extractProductList {
    my($foldList,$distinctReadsHash,$parameters) = @_;
    my %productHash;
    foreach my $foldInfo (@{$foldList}) {
	my($name,$location,$mfe,$seq,$fold) = @{$foldInfo};
	my $leftCenter = getMaxLeftCenter($fold);
	my $rightCenter = getMaxRightCenter($fold);
	my $basePairs = getBasePairs($leftCenter,$rightCenter,
				     $fold,$parameters);
	my $productInfo = extractProducts($leftCenter,$rightCenter,$fold,$basePairs,$location,
					  $distinctReadsHash->{$name},$parameters);
	$productHash{$name} = $productInfo;
    }
    return \%productHash;
}

sub extractProducts {
    my($leftCenter,$rightCenter,$fold,$basePairs,$location,$distinctReads,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $productInfo = getProductInfo($leftCenter,$rightCenter,$fold,$basePairs,
				     $start,$stop,$strand,
				     $distinctReads,$parameters);
    $productInfo = rebuildProducts($leftCenter,$rightCenter,$fold,$basePairs,
				    $productInfo,$parameters);    
    return $productInfo;
}

sub getMirCount {
    my($productInfo) = @_;
    my $mirCount = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};	
	if($newType eq "miR") {
	    $mirCount += $prodCount;
	}
    }
    return $mirCount;
}


sub computeMaxSameShift {
    my($leftCenter,$rightCenter,$fold,$basePairs,$location,$distinctReads,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $newParameters = {};
    foreach my $key (keys %{$parameters}) {
	$newParameters->{$key} = $parameters->{$key};
    }
    $newParameters->{minDist} = $parameters->{minShift};
    my $totalReads = getTotalReads($distinctReads);
    #print "total reads: $totalReads\n";
    my $productInfo = getProductInfo($leftCenter,$rightCenter,$fold,$basePairs,
				     $start,$stop,$strand,
				     $distinctReads,$newParameters);
    #print "product info:\n";
    foreach my $products (@{$productInfo}) {
	# note, these products are not "rebuilt" with rebuildProducts
	my($side,$productList,$productCount,$maxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$products};
	#print "$side $relStart $relStop : $productCount $maxProductCount\n";
    }
    $productInfo = rebuildProducts($leftCenter,$rightCenter,$fold,$basePairs,
				   $productInfo,$newParameters);
    return getSameArmShift($productInfo,$totalReads,$newParameters);
}

sub computeMaxBothShift {
    my($leftCenter,$rightCenter,$fold,$basePairs,$location,$distinctReads,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $newParameters = {};
    foreach my $key (keys %{$parameters}) {
	$newParameters->{$key} = $parameters->{$key};
    }
    $newParameters->{minDist} = $parameters->{minShift};
    my $totalReads = getTotalReads($distinctReads);
    #print "total reads: $totalReads\n";
    my $productInfo = getProductInfo($leftCenter,$rightCenter,$fold,$basePairs,
				     $start,$stop,$strand,
				     $distinctReads,$newParameters);
    #print "product info:\n";
    foreach my $products (@{$productInfo}) {
	# note, these products are not "rebuilt" with rebuildProducts
	my($side,$productList,$productCount,$maxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$products};
	#print "$side $relStart $relStop : $productCount $maxProductCount\n";
    }
    $productInfo = rebuildProducts($leftCenter,$rightCenter,$fold,$basePairs,
				   $productInfo,$newParameters);
    if(bothArmProducts($basePairs,$productInfo,$newParameters)) {
	#print "possible both arm shift\n";
	return getBothArmShift($basePairs,$productInfo,$totalReads,$newParameters,$location);
    }
    return 0;
}

sub bothArmProducts {
    my($basePairs,$productInfo,$parameters) = @_;
    my %sideHash;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};	
	if($side eq "3p") {
	    if(overlapsBasePairs($basePairs,$relStart,$relStop)) {
		$sideHash{$side}++;
	    }
	} else {
	    $sideHash{$side}++;
	}
    }
    if(($sideHash{"5p"})&&($sideHash{"3p"})) {
	return 1;
    }
    return 0;
}

sub getBothArmShift {
    my($basePairs,$productInfo,$totalReads,$parameters,$location) = @_;
    my $minReadFraction = 0.01;
    my $countThreshold = $totalReads*$minReadFraction;
    my $minOverlap = $parameters->{minOverlap} or die "FAIL: no minOverlap loaded.\n";
    my $maxShift = 0;
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$prodCount1,$maxProdCount1,
	   $relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	if($prodCount1 >= $countThreshold) {	
	    if($side1 eq "5p") {
		for(my $j=0;$j<@{$productInfo};$j++) {
		    if($i != $j) {
			my($side2,$newType2,$prodList2,$prodCount2,$maxProdCount2,
			   $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			if($prodCount2 >= $countThreshold) {
			    if($side2 eq "3p") {				
				if(overlapsBasePairs($basePairs,$relStart2,$relStop2)) {
				    my($relStart2bp,$relStop2bp) = getOuterBasePairsFivePrime($basePairs,
											      $relStart2,
											      $relStop2,
											      $location);
				    my $overlap = getOverlap($relStart1,$relStop1,
							     $relStart2bp,$relStop2bp);
				    my $shift = getShift($relStart1,$relStop1,
							 $relStart2bp,$relStop2bp);
				    if($overlap > $minOverlap) {
					if(abs($shift) > abs($maxShift)) {				
					    $maxShift = $shift;
					}
				    }
				}	
			    }
			}
		    }
		}
	    }
	}
    }
    return $maxShift;
}

sub overlapsBasePairs {
    my($basePairs,$relStart,$relStop) = @_;
    my @pairs;
    foreach my $pair (@{$basePairs}) {
	my($bp1,$bp2) = @{$pair};
	if(($relStart <= $bp2)&&($bp2 <= $relStop)) {
	    push(@pairs,[$bp1,$bp2]);
	}
    }
    if(@pairs > 1) {
	return 1;
    } else {
	return 0;
    }
}

sub computeBasePairDensity {
    my($basePairs,$relStart,$relStop) = @_;
    my $basePairCount = 0;
    foreach my $pair (@{$basePairs}) {
	my($bp1,$bp2) = @{$pair};
	if(($relStart <= $bp2)&&($bp2 <= $relStop)||
	   ($relStart <= $bp1)&&($bp1 <= $relStop)) {	    
	    $basePairCount++;
	}
    }
    return sprintf("%.2f",$basePairCount / ($relStop - $relStart + 1));
}

sub getOuterBasePairs {
    # NOTE: This assumes that relStart and relStop are on the 3p arm!
    my($basePairs,$relStart,$relStop,$location) = @_;
    my @pairs;
    foreach my $pair (@{$basePairs}) {
	my($bp1,$bp2) = @{$pair};
	# here assumes bp2 is basepair on 3p arm
	if(($relStart <= $bp2)&&($bp2 <= $relStop)) {
	    push(@pairs,[$bp1,$bp2]);
	}
    }
    if(@pairs > 1) {
	# sort by bp position on 5' arm
	@pairs = sort {$a->[0] <=> $b->[0]} @pairs;
	my @runList;
	my $runStart = 0;
	my $runStop = 0;
	my $runStartBp = 0;
	my $runStopBp = 0;
	my $INRUN = 0;
	for(my $i=0;$i<@pairs-1;$i++) {
	    my($bp1,$bp2) = @{$pairs[$i]};
	    my($nbp1,$nbp2) = @{$pairs[$i+1]};	    
	    # if base pairs on both arms are adjacent nucleotides.
	    if(($bp1+1==$nbp1)&&($bp2-1==$nbp2)) {
		# current base pairs are neighbors.
		if($INRUN) {
		    # in a run of neighboring base pairs, extend the run
		    $runStop = $nbp1;
		    $runStopBp = $nbp2;
		} else {
		    # not in a run, start a new run
		    $runStart = $bp1;
		    $runStop = $nbp1;
		    $runStartBp = $bp2;
		    $runStopBp = $nbp2;
		    $INRUN = 1;
		}
	    } else {
		if($INRUN) {
		    push(@runList,[$runStart,$runStop,$runStartBp,$runStopBp]);
		}
		$INRUN = 0;
	    }
	}
	if($INRUN) {
	    push(@runList,[$runStart,$runStop,$runStartBp,$runStopBp]);
	}
	my $maxLength = 0;
	my $leftPair;
	my $rightPair;
	# here runStart is the most 5' basepair that pairs with a base that overlaps
	# the relStart and relStop.
	foreach my $run (@runList) {
	    my($runStart,$runStop,$runStartBp,$runStopBp) = @{$run};
	    if($runStop - $runStart > $maxLength) {
		$maxLength = $runStop - $runStart;
		$leftPair = [$runStart,$runStartBp];
		$rightPair = [$runStop,$runStopBp];
	    }
	}
	if($maxLength) {
	    unless($leftPair) {
		my $L = "L" x ($relStart+1);
		my $R = "R" x ($relStop+1);
		die "could not get leftPair:\n$L\n$R\nat $location\n";
	    }
	    # lbp1 is the most 5' base overlapping with relStart,relStop
	    my($lbp1,$lbp2) = @{$leftPair};
	    my $deltaL = $relStop - $lbp2;
	    unless($rightPair) {
		my $L = "L" x ($relStart+1);
		my $R = "R" x ($relStop+1);
		die "could not get rightPair:\n$L\n$R\nat $location\n";
	    }
	    # lbp1 is the most 3' base overlapping with relStart,relStop
	    my($rbp1,$rbp2) = @{$rightPair};
	    my $deltaR = $rbp2 - $relStart;
	    if(($deltaL < 0)||($deltaR < 0)) {
		die "left: $lbp1 $lbp2\nright:$rbp1 $rbp2\nwindow: $relStart $relStop\n";
	    }
	    #print "longest run: $lbp1,$lbp2 .. $rbp1,$rbp2\n";
	    return ($lbp1-$deltaL,$rbp1+$deltaR);
	}
    }
    return (0,0);
}

sub getOuterBasePairsFivePrime {
    # NOTE: This assumes that relStart and relStop are on the 3p arm!
    my($basePairs,$relStart,$relStop,$location) = @_;
    my @pairs;
    foreach my $pair (@{$basePairs}) {
	my($bp1,$bp2) = @{$pair};
	# here assumes bp2 is basepair on 3p arm
	if(($relStart <= $bp2)&&($bp2 <= $relStop)) {
	    push(@pairs,[$bp1,$bp2]);
	}
    }
    # sort by bp position on 5' arm
    if(@pairs > 1) {
	@pairs = sort {$a->[0] <=> $b->[0]} @pairs;
	# lbp1 is the most 5' base overlapping with relStart,relStop
	my $leftPair = shift(@pairs);
	my $rightPair = pop(@pairs);
	my($lbp1,$lbp2) = @{$leftPair};
	my $deltaL = $relStop - $lbp2;
	unless($rightPair) {
	    my $L = "L" x ($relStart+1);
	    my $R = "R" x ($relStop+1);
	    die "could not get rightPair:\n$L\n$R\nat $location\n";
	}
	# lbp1 is the most 3' base overlapping with relStart,relStop
	my($rbp1,$rbp2) = @{$rightPair};
	my $deltaR = $rbp2 - $relStart;
	if(($deltaL < 0)||($deltaR < 0)) {
	    die "left: $lbp1 $lbp2\nright:$rbp1 $rbp2\nwindow: $relStart $relStop\n";
	}
	#print "longest run: $lbp1,$lbp2 .. $rbp1,$rbp2\n";
	return ($lbp1-$deltaL,$rbp1+$deltaR);
    }
    return (0,0);
}

sub getSameArmShift {
    my($productInfo,$totalReads,$parameters) = @_;
    my $minOverlap = $parameters->{minOverlap} or die "FAIL: no minOverlap loaded.\n";
    my $minReadFraction = 0.01;
    my $countThreshold = $totalReads*$minReadFraction;
    my $maxShift = 0;
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$prodCount1,$maxProdCount1,
	   $relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	#print "$i: $side1 $newType1 $relStart1 $relStop1 $prodCount1 vs $countThreshold\n";
	if($prodCount1 >= $countThreshold) {	    
	    for(my $j=0;$j<@{$productInfo};$j++) {
		if($i != $j) {
		    my($side2,$newType2,$prodList2,$prodCount2,$maxProdCount2,
		       $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
		    if($prodCount2 >= $countThreshold) {
			my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2);
			my $shift = getShift($relStart1,$relStop1,$relStart2,$relStop2);
			if($overlap > $minOverlap) {
			    if($shift > $maxShift) {
				$maxShift = $shift;				
			    }
			}
		    }
		}
	    }
	}
    }
    return $maxShift;
}

sub lociContained {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    if(($relStart1 <= $relStart2)&&($relStop2 <= $relStop1)) {
	return 1;
    } elsif(($relStart2 <= $relStart1)&&($relStop1 <= $relStop2)) {
	return 1;
    }
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

sub getShift {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    return $relStart2-$relStart1;
}

sub doubleStrandedProducts {
    my($productInfo,$revProductInfo,$parameters) = @_;
    my $maxReverse = $parameters->{maxReverse};
    my $readCount = getTotalProductReads($productInfo);
    my $revReadCount = getTotalProductReads($revProductInfo);
    my $totalReads = $readCount + $revReadCount;
    my $fraction = $totalReads ? $revReadCount/$totalReads : 0;
    #print "doubleStrandedProducts $totalReads $readCount $revReadCount $fraction\n";
    if($fraction <= $maxReverse) {
	return 0;
    } else {
	return 1;
    }
}

sub overlappingProducts {
    my($productInfo,$revProductInfo,$parameters) = @_;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	my $start = $relStart;
	my $stop = $relStop;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aProdCount,$aMaxProdCount,
	       $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
	    my $aStart = $aRelStart;
	    my $aStop = $aRelStop;
	    # if the products overlap..
	    if((($start <= $aStart)&&($aStart <= $stop))||
	       (($start <= $aStop)&&($aStop <= $stop))||
	       (($aStart <= $start)&&($start <= $aStop))||
	       (($aStart <= $stop)&&($stop <= $aStop))) {
		return 1;
	    }
	}
    }
    # if made it here, then no overlaps.
    return 0;
}

sub getReverseProductDisplacement {
    my($productInfo,$revProductInfo) = @_;
    my $minDispInit = 100;
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    my $totalReads = getTotalProductReads($productInfo);
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	my $start = $relStart;
	my $stop = $relStop;
	my $minDisp = $minDispInit;
	my $minI2;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aProdCount,$aMaxProdCount,
	       $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
	    if($aProdCount >= $minFrac*$totalReads) {
		my $aStart = $aRelStart;
		my $aStop = $aRelStop;
		# if the products overlap..
		if((($start <= $aStart)&&($aStart <= $stop))||
		   (($start <= $aStop)&&($aStop <= $stop))||
		   (($aStart <= $start)&&($start <= $aStop))||
		   (($aStart <= $stop)&&($stop <= $aStop))) {
		    my $disp = abs($aStart - $start);
		    unless($USED{$i2}) {
			if($disp < $minDisp) {
			    $minDisp = $disp;
			    $minI2 = $i2;
			}
		    }
		}
	    }
	}
	if($minDisp < $minDispInit) {
	    $sum += $minDisp;
	    $total++;
	    $USED{$minI2}++;	    
	}
    }
    return($sum,$total);
}

sub goodCandidate {
    my($id,$location,$tRNAs,$parameters) = @_;
    if($tRNAs->{$id}) {
	open(FRC,">>$parameters->{filteredCandidatesFile}");
	print FRC "$id\t$location\trejected: tRNA overlap\n";
	close(FRC);
	return 0;
    }
    return 1;
}

sub plausibleReads {
    my($productInfo,$revProductInfo,$parameters) = @_;
    if(doubleStrandedProducts($productInfo,$revProductInfo,$parameters)) {
	if(overlappingProducts($productInfo,$revProductInfo,$parameters)) {
	    return (1,"");
	} else {
	    return (0,"rejected: has abundant reads on each strand, but none overlap.");
	}
    } else {
	return (1,"");
    }
}

sub computeProductBasePairing {
    my($leftCenter,$rightCenter,$productInfo,$basePairs,$parameters) = @_;
    my $minBasePairDensity = 1.0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	if($newType eq "miR") {	    
	    my $basePairDensity = computeBasePairDensity($basePairs,$relStart,$relStop);
	    if($basePairDensity < $minBasePairDensity) {
		$minBasePairDensity = $basePairDensity;
	    }
	}
    }
    return $minBasePairDensity;    
}

sub goodProducts {
    my($productInfo,$parameters) = @_;
    my $maxFivePrimeHet = $parameters->{maxFivePrimeHet} or die "FAIL: no maxFivePrimeHet loaded.\n";
    my $minLocusCount = $parameters->{minLocusCount} or die "FAIL: no minLocusCount loaded.\n";
    my $readCount = getTotalProductReads($productInfo);
    #print "readCount = $readCount\n";
    if($readCount < $minLocusCount) {
	return 0;
    }
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount) = @{$product};
	my $fivePrimeHet = computeFivePrimeHet($prodList);
	#print "prodCount: $prodCount, 5' het: $fivePrimeHet\n";
	if(($fivePrimeHet < $maxFivePrimeHet)&&($prodCount > 1)) {
	    return 1;
	}
    }
    return 0;
}

sub computeMaxProdFivePrimeHet {
    my($productInfo,$parameters) = @_;
    my $maxCount = 0;
    my $maxProdFivePrimeHet = 0;
    foreach my $product (@{$productInfo}) {
        my($side,$newType,$prodList,$prodCount,$maxProdCount) = @{$product};
	if($newType eq "miR") {
	    my $fivePrimeHet = computeFivePrimeHet($prodList);
	    if($prodCount > $maxCount) {
		$maxProdFivePrimeHet = $fivePrimeHet;
		$maxCount = $maxProdCount;
	    }
	}
    }
    return $maxProdFivePrimeHet;
}

sub computeAverageFivePrimeHet {
    my($productInfo,$parameters) = @_;
    my $sum = 0.0;
    my $total = 0;
    foreach my $product (@{$productInfo}) {
        my($side,$newType,$prodList,$prodCount,$maxProdCount) = @{$product};
        my $fivePrimeHet = computeFivePrimeHet($prodList);
	$sum += $fivePrimeHet;
	$total++;
    }
    return $total ? $sum/$total : 0.0;
}

sub computeFivePrimeHet {
    my($productList) = @_;
    my $FIRST=1;
    my $fivePrimeMaxPos;
    my $fivePrimeMaxCount = 0;
    my $fivePrimeTotalCount = 0;
    my %startCount;
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$libraryCounts) = @{$read};
	$startCount{$relStart} += $count;
    }
    my @starts = sort {$startCount{$b} <=> $startCount{$a}} keys(%startCount);
    my $topStart = shift(@starts);
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$libraryCounts) = @{$read};
	if($relStart == $topStart) {
	    $fivePrimeMaxCount += $count;
	}
        $fivePrimeTotalCount += $count;
    }
    my $fivePrimeHet = ($fivePrimeTotalCount-$fivePrimeMaxCount)/$fivePrimeTotalCount;
    return $fivePrimeHet;
}

sub rebuildProducts {
    my($leftCenter,$rightCenter,$fold,$basePairs,$productInfo,$parameters) = @_;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
    my %usedMirSide;
    my %usedMorSide;
    my %newTypes;
    my @newProductInfo1;
    my @newProductInfo2;
    # sort by distance to loop, ascending order.
    my $totalReads = 0;
    foreach my $products (@{$productInfo}) {
	my($side,$productList,$productCount,$maxProductCount) = @{$products};
	$totalReads += $productCount;
    }   
    my @sortedProductInfo = sort {abs($a->[6]) <=> abs($b->[6])} @{$productInfo};
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$productCount,$maxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$products};
	my $newType;
	#print "$side $relStart $relStop $offset\n";
	if($side eq "loop") {
	    $newType = "loop";
	} elsif($side eq "split") {
	    $newType = "split";
	} elsif($side eq "out") {
	    $newType = "out";
	} elsif(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if(onHairpinArm($leftCenter,$rightCenter,$leftEnd,$rightEnd,
			    $relStart,$relStop,$parameters)) {		
		if($usedMirSide{$side}) {
		    if($usedMorSide{$side}) {
			$newType = "out";
		    } else {
			$newType = "moR";
			$usedMorSide{$side}++;
		    }
		} elsif($usedMirSide{otherSide($side)}) {
		    if(overlapsMir($leftCenter,$rightCenter,$basePairs,$side,$relStart,$relStop,\@newProductInfo1,$parameters)) {
			$newType = "miR";
			$usedMirSide{$side}++;
		    } else {
			if($usedMorSide{$side}) {
			    $newType = "out";
			} else {
			    $newType = "moR";
			    $usedMorSide{$side}++;
			}
		    }
		} else {
		    if($productCount > 0.05*$totalReads) {
			$newType = "miR";
			$usedMirSide{$side}++;
		    } else {
			$newType = "loop";
		    }
		}
	    } else {
		# this should be redundant, but here just in case.
		$newType = "loop";
	    }
	} elsif(closeToHairpin($leftCenter,$rightCenter,$relStart,$relStop,$parameters)) {
	    if($usedMirSide{$side}) {
		if($usedMorSide{$side}) {
		    $newType = "out";
		} else {
		    $newType = "moR";
		    $usedMorSide{$side}++
		}
	    } elsif($usedMirSide{otherSide($side)}) {
		if(overlapsMir($leftCenter,$rightCenter,$basePairs,
			       $side,$relStart,$relStop,\@newProductInfo1,$parameters)) {
		    $newType = "miR";
		    $usedMirSide{$side}++;
		} else {
		    if($usedMorSide{$side}) {
			$newType = "out";
		    } else {
			$newType = "moR";
			$usedMorSide{$side}++;
		    }
		}
	    } else {
		$newType = "out";
	    }
	} else {
	    $newType = "out";
	}
	$newTypes{$relStart."x".$relStop} = $newType;
	push(@newProductInfo1,[$side,$newType,$productList,
			      $productCount,$maxProductCount,
			      $relStart,$relStop,$offset,$gStart]);
    }
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$productCount,$maxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$products};
	if(($side eq "5p")||($side eq "3p")) {
	    if($newTypes{$relStart."x".$relStop} eq "loop") {
		if(overlapsMir($leftCenter,$rightCenter,$basePairs,$side,
			       $relStart,$relStop,\@newProductInfo1,$parameters)) {
		    $newTypes{$relStart."x".$relStop} = "miR";
		}
	    }
	}
	push(@newProductInfo2,[$side,$newTypes{$relStart."x".$relStop},$productList,
			       $productCount,$maxProductCount,
			       $relStart,$relStop,$offset,$gStart]);
    }
    return \@newProductInfo2;
}

sub otherSide {
    my $side = shift;
    if($side eq "5p") {
	return "3p";
    } elsif($side eq "3p") {
	return "5p";
    } else {
	die "unexpected side in otherSide: $side\n";
    }
}

sub overlapsMir {
    my($leftCenter,$rightCenter,$basePairs,$side,$relStart,$relStop,$newProductInfo,$parameters) = @_;
    my $minShift = $parameters->{minShift} or die "FAIL: no minShift loaded.\n";
    foreach my $oProduct (@{$newProductInfo}) {
	my($oSide,$oNewType,$oProductlist,$oProductCount,$oMaxProdCount,$oRelStart,$oRelStop) = @{$oProduct};
	if($oNewType eq "miR") {
	    if(($side eq "5p")&&($oSide eq "3p")) {
		my($oRelStartBp,$oRelStopBp) = getOuterBasePairs($basePairs,$oRelStart,$oRelStop);		
		my $overlap = getOverlap($relStart,$relStop,$oRelStartBp,$oRelStopBp);
		my $shift = getShift($relStart,$relStop,$oRelStartBp,$oRelStopBp);
		#print "1:$side ($relStart,$relStop) ($oRelStartBp,$oRelStopBp) $shift\n";
		if($overlap > 0) { 
		    if(abs($shift) <= $minShift) {
			return 1;
		    }
		}
	    } elsif(($side eq "3p")&&($oSide eq "5p")) {
		my($relStartBp,$relStopBp) = getOuterBasePairs($basePairs,$relStart,$relStop);
		my $shift = getShift($oRelStart,$oRelStop,$relStartBp,$relStopBp);
		my $overlap = getOverlap($oRelStart,$oRelStop,$relStartBp,$relStopBp);
		#print "2:$side ($relStartBp,$relStopBp) ($oRelStart,$oRelStop) $shift\n";
		if($overlap > 0) {
		    if(abs($shift) <= $minShift) {
			return 1;
		    }
		}
	    }
	}
    }
    return 0;
}

sub withinHairpin {
    my($leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;
    my $inBuffer = $parameters->{inHairpinBuffer};
    my $outBuffer = $parameters->{outHairpinBuffer};   
    if(($leftEnd-$outBuffer<=$relStart)&&
       ($relStop<=$rightEnd+$outBuffer)) {	    
	return 1;
    }
    return 0;
}

sub onHairpinArm {
    my($leftCenter,$rightCenter,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;
    my $inBuffer = $parameters->{inHairpinBuffer};
    my $outBuffer = $parameters->{outHairpinBuffer};   
    if((($leftEnd-$outBuffer<=$relStart)&&($relStop<=$leftCenter+$inBuffer))||
       (($rightCenter-$inBuffer<=$relStart)&&($relStop<=$rightEnd+$outBuffer))) {	    
	return 1;
    }
    return 0;
}

sub closeToHairpin {
    my($leftCenter,$rightCenter,$relStart,$relStop,$parameters) = @_;
    my $hairpinRange = $parameters->{hairpinRange} or die "no hairpinRange loaded\n";
    if(($leftCenter-$hairpinRange <= $relStart)&&
       ($relStop <= $rightCenter+$hairpinRange)) {
	return 1;
    }
    return 0;
}

sub doubleStranded {
    my($start,$stop,$strand,$chromDataSet,$hitCount,$parameters) = @_;
    my $maxReverse = $parameters->{maxReverse};
    my $senseReadCount = readCount($start,$stop,$chromDataSet->{$strand},$hitCount);
    my $revStrand = revStrand($strand);
    my $reverseReadCount = readCount($start,$stop,$chromDataSet->{$revStrand},$hitCount);
    my $readCount = $senseReadCount + $reverseReadCount;
    # if reverse is less than 5%, it is singleStranded
    my $fraction = $readCount ? $reverseReadCount/$readCount : 0;
    #print "doubleStranded $readCount $senseReadCount $reverseReadCount $fraction\n";
    if($fraction <= $maxReverse) {
	return 0;
    } else {
	return 1;
    }
}

sub computeMaxProdHitCount {
    my($productInfo,$location,$dataSet,$hitCount,$parameters) = @_;
    my $minDist = $parameters->{minDist};
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $maxCount = 0;
    my $maxProdHitCount = 0;
    foreach my $product (@{$productInfo}) {
        my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	if($newType eq "miR") {
	    my $sum = 0;
	    my $total = 0;
	    foreach my $read (@{$dataSet->{$chrom}{$strand}}) {
		# for each read...
		my($rStart,$rStop,$sample,$rID) = @{$read};
		# if it is associated with a product...
		if(abs($gStart-$rStart) < $minDist) {
		    # compute the average number of hits to the genome
		    $sum += $hitCount->{$sample.$rID};
		    $total++;
		    #print "incrementing: ( $rStart $rStop ) $rID $hitCount->{$sample.$rID}\n";
		    #print "--> $sum $total\n";
		}
	    }
	    my $avg = $total ? $sum/$total : 0.0;
	    #print " = $avg\n";
	    if($prodCount > $maxCount) {
		$maxCount = $prodCount;
		$maxProdHitCount = $avg;
		
	    }
	}
    }
    return $maxProdHitCount;
}

sub computeMaxProdHitCount2 {
    my($productInfo,$distinctReads,$parameters) = @_;
    my $minDist = $parameters->{minDist};
    my $maxCount = 0;
    my $maxProdHitCount = 0;
    foreach my $product (@{$productInfo}) {
        my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	if($newType eq "miR") {
	    my $sum = 0;
	    my $total = 0;
	    foreach my $dRead (@{$distinctReads}) {
		# for each read...
		my($dStart,$dStop,$dTotal,$libraryCounts,$dTotalHitCount,$dReadTotal) = @{$dRead};
		# if it is associated with a product...
		if(abs($gStart-$dStart) < $minDist) {
		    # compute the average number of hits to the genome
		    $sum += $dTotalHitCount;
		    $total += $dReadTotal;
		}
	    }	
	    my $avg = $total ? $sum/$total : 0.0;
	    if($prodCount > $maxCount) {
		$maxCount = $prodCount;
		$maxProdHitCount = $avg;
		
	    }
	}
    }
    return $maxProdHitCount;
}

sub computeUniqueReadFraction {
    my($id,$location,$dataSet,$hitCount) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $uniqueReadCount = uniqueReadCount($start,$stop,$dataSet->{$chrom}{$strand},$hitCount);
    my $readCount = readCount($start,$stop,$dataSet->{$chrom}{$strand},$hitCount);
    return $readCount ? $uniqueReadCount/$readCount : 0;
}

sub computeUniqueReadFraction2 {
    my($distinctReads) = @_;
    my $uniqueReadCount = 0;
    my $readCount = 0;
    foreach my $dRead (@{$distinctReads}) {
	my($dStart,$dStop,$dTotal,$libraryCounts,$totalHitCount,$dReadTotal) = @{$dRead};
	if($totalHitCount == $dReadTotal) {
	    $uniqueReadCount += $dTotal;
	}
	$readCount += $dTotal;
    }
    return $readCount ? $uniqueReadCount/$readCount : 0;
}

sub uniqueReadCount {
    my($start,$stop,$chromStrandDataSet,$hitCount) = @_;
    my $uniqueReadCount = 0;
    foreach my $read (@{$chromStrandDataSet}) {
        my($rStart,$rStop,$sample,$rID) = @{$read};
        if(($start <= $rStart)&&($rStop <= $stop)) {
            if($hitCount->{$sample.$rID} == 1) {
		$uniqueReadCount++;
	    }
        }
    }
    return $uniqueReadCount;
}

sub readCount {
    my($start,$stop,$chromStrandDataSet,$hitCount) = @_;
    my $readCount = 0;
    foreach my $read (@{$chromStrandDataSet}) {
        my($rStart,$rStop,$sample,$rID) = @{$read};
        if(($start <= $rStart)&&($rStop <= $stop)) {
            $readCount += 1/($hitCount->{$sample.$rID});
        }
    }
    return $readCount;
}

sub revStrand {
    my $strand = shift;
    if($strand eq "-") {
	return "+";
    }
    return "-";
}

sub getTotalReads {
    my($distinctReads) = @_;
    my $totalReads = 0;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$libraryCounts) = @{$dRead};
	$totalReads += $total;
    }
    return $totalReads;
}

sub getProductInfo {
    my($leftCenter,$rightCenter,$fold,$basePairs,
       $start,$stop,$strand,$distinctReads,$parameters) = @_;
    # first build a product hash:
    my $PRODUCTCOUNT = 1;
    my %productHash;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$libraryCounts) = @{$dRead};
	my $relStart = $dStart - $start;
	my $relStop = $dStop - $start;
	# relStart is the position relative to the 5' end of readRegion. numbers start from zero.
	if($strand eq "-") {
	    $relStart = $stop - $dStop;
	    $relStop = $stop - $dStart;
	}
	# offset is relative to the leftMiddlePos (last 5' base pair in hairpin)
	my $offset = $relStop - $leftCenter;	
	if(($offset > 0)&&($relStop < $rightCenter)) {
	    $offset = 0;
	} elsif($offset > 0) {
	    $offset = $relStart - $rightCenter;
	}
	my $parsedRead = [$relStart,$relStop,$offset,$dStart,
			  $total,$libraryCounts];
	if(my $id = overlapCurrent(\%productHash,$parsedRead,$parameters)) {
	    push(@{$productHash{$id}},$parsedRead);	    
	} else {
	    push(@{$productHash{$PRODUCTCOUNT}},$parsedRead);
	    $PRODUCTCOUNT++;
	}
    }
    my @productList;
    foreach my $id (keys %productHash) {
	my $total=0;
	my $productStart;
	my $productStop;
	my $FIRST=1;
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$offset,$dStart,$readCount,$libraryCounts) = @{$storedRead};
	    $total += $readCount;
	    if($FIRST) {
		$FIRST=0;  
		$productStart = $relStart;
		$productStop = $relStop;
	    }
	}
	push(@productList,[$id,$total,$productStart,$productStop]);
    }
    # now sort the products, by total, highest to lowset.                    
    my @sortedProductList = sort {$b->[1] <=> $a->[1]} @productList;
    my @productInfo;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
    foreach my $idInfo (@sortedProductList) {
	my($id,$totalCount,$productStart,$productStop) = @{$idInfo};
	# get the relative position that is most abundant:
	my($maxRelStart,$maxRelStop,$offset,$gStart) = getMaxProductPos($productHash{$id});
	my($maxProductCount,$productCount) = getProductCounts($productHash{$id});
	my $side = getProductSide($leftCenter,$rightCenter,$leftEnd,$rightEnd,
				  $maxRelStart,$maxRelStop,$parameters);
	push(@productInfo,[$side,$productHash{$id},$productCount,$maxProductCount,
			   $maxRelStart,$maxRelStop,$offset,$gStart]);
    }
    return \@productInfo;
}

sub getOverlappingReadCount {
    my($chrom,$start,$stop,$strand,$dataSet,$hitCount) = @_;
    my $readCount = 0;
    foreach my $read (@{$dataSet->{$chrom}{$strand}}) {
	my($rStart,$rStop,$sample,$rID) = @{$read};
	if(($start <= $rStart)&&($rStop <= $stop)) {
	    $readCount += 1/($hitCount->{$sample.$rID});	    
	}	
    }
    return $readCount;
}

sub extractDistinctReadInfo {
    my($start,$stop,$chromStrandDistinctReads) = @_;
    my @distinctReads;
    foreach my $dRead (@{$chromStrandDistinctReads}) {
	my($dStart,$dStop,$total,$libraryTotals,$totalHitCount,$dReadTotal) = @{$dRead};
	if(($start<=$dStart)&&($dStop<=$stop)) {
	    push(@distinctReads,[$dStart,$dStop,$total,$libraryTotals,$totalHitCount,$dReadTotal]);
	}
    }
    @distinctReads = sort {$b->[2] <=> $a->[2]} @distinctReads;
    return \@distinctReads;
}

sub getDistinctReadInfo {
    my($start,$stop,$chromStrandDataSet,$sampleList,$hitCount) = @_;
    my %distinctReadCount;
    my %libraryCounts;
    foreach my $read (@{$chromStrandDataSet}) {
	my($rStart,$rStop,$sample,$rID) = @{$read};
	if(($start-1 <= $rStart)&&($rStop <= $stop+1)) {
	    $distinctReadCount{$rStart."x".$rStop} += 1/($hitCount->{$sample.$rID});
	    $libraryCounts{$rStart."x".$rStop}{$sample} += 1/($hitCount->{$sample.$rID});
	}
    }
    my @distinctReads;
    foreach my $key (keys %distinctReadCount) {
	my($drStart,$drStop) = split(/x/,$key);
	foreach my $sample (@{$sampleList}) {
	    unless($libraryCounts{$key}{$sample}) {
		$libraryCounts{$key}{$sample} = 0;
	    }
	}
	push(@distinctReads,[$drStart,$drStop,$distinctReadCount{$key},$libraryCounts{$key}]);
    }
    # sort distinct reads by abundance.  This is important in getProductInfo
    @distinctReads = sort {$b->[2] <=> $a->[2]} @distinctReads;
    return \@distinctReads;
}

sub overlapCurrent {
    my($productHash,$parsedRead,$parameters)=@_;
    my $minDist = $parameters->{minDist};
    my($readRelStart,$readRelStop,$offset,$gStart,$count,$libraryCounts) = @{$parsedRead};
    # sort by the first entries total count ( the most abundant entry by construction, since
    # distinctReads is already sorted by abundance. ) highest to lowest...
    my @sortedIdList = sort {$productHash->{$b}[0][4] <=> $productHash->{$a}[0][4]} keys(%{$productHash});
    my @idList;
    foreach my $id (@sortedIdList) {
	my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount)=@{$productHash->{$id}[0]};
	#print "comparing($id): ($readRelStart,$readRelStop) to ($sRelStart,$sRelStop) : ",abs($sRelStart-$readRelStart),"\n";
	if(abs($sRelStart-$readRelStart) <= $minDist) {
	    push(@idList,$id);
	    #print "ACCPETED -> $id\n";
	}
    }
    if(scalar(@idList) == 0) {
	# no overlaps, return FALSE.
        return "";
    } else {
	# more than one overlap, return the most abundant, first one by construction..
        return shift(@idList);
    }
}

sub getProductSide {
    my($leftCenter,$rightCenter,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;    
    my $inBuffer = $parameters->{inHairpinBuffer} or die "no inHairpinBuffer loaded\n";
    my $outBuffer = $parameters->{outHairpinBuffer} or die "no outHairpinBuffer loaded\n";   
    my $hairpinRange = $parameters->{hairpinRange} or die "no hairpinRange loaded";
    if(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	if(onHairpinArm($leftCenter,$rightCenter,$leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if($relStop <= $leftCenter+$inBuffer) {
		return "5p";
	    } elsif($rightCenter-$inBuffer<=$relStart) {	    
		return "3p";	
	    } else {
		die "unknow situation reached in getProductSide!!! ouch.\n";
	    }
	} elsif(($relStart<=$leftCenter+1)&&($rightCenter-1<=$relStop)) {
	    return "loop";
	} elsif(($leftCenter<=$relStart)&&($relStop<=$rightCenter)) {
	    return "loop";
	} else {
	    return "split";
	}
    } elsif(($leftCenter-$hairpinRange <= $relStart)&&
	    ($relStop<=$rightCenter+$hairpinRange)) {
	if(($leftCenter-$hairpinRange <= $relStart)&&
	   ($relStop <= $leftCenter+$inBuffer)) {
	    return "5p";
	} elsif(($rightCenter-$inBuffer<=$relStart)&&
		($relStop<=$rightCenter+$hairpinRange)) {	    
	    return "3p";
	}	
    } else {
	return "out";
    }
    return "out";
}

sub getMaxProductPos {
    my $productList = shift;
    my($maxRelStart,$maxRelStop,$maxOffset,$maxGStart);
    my $FIRST = 1;
    # list is assumed to be sorted by abundance. must be sorted in getProductInfo()
    foreach my $read (@{$productList}) {
        my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount)=@{$read};
        if($FIRST) {
            $FIRST=0;
            $maxRelStart = $sRelStart;
	    $maxRelStop = $sRelStop;
	    $maxOffset = $sOffset;
	    $maxGStart = $sGStart;
	}
    }
    return($maxRelStart,$maxRelStop,$maxOffset,$maxGStart);
}

sub getTotalProductReads {
    my $productInfo = shift;
    my $totalReads = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	$totalReads += $prodCount;
    }
    return $totalReads;
}

sub getProductCounts {
    my $productList = shift;
    my $productCounts = 0;
    my $maxProductCount = 0;
    my $FIRST = 1;
    foreach my $read (@{$productList}) {
        my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount)=@{$read};
        $productCounts += $sCount;
        if($FIRST) {
            $FIRST=0;
            $maxProductCount = $sCount;
        }
    }
    return($maxProductCount,$productCounts);
}

sub getTotalLibraryCounts {
    my $productList = shift;
    my %totalLibraryCounts;
    foreach my $read (@{$productList}) {
	my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount,$libraryCounts)=@{$read};
	foreach my $sample (keys %{$libraryCounts}) {
	    $totalLibraryCounts{$sample} += $libraryCounts->{$sample};
	}
    }
    return \%totalLibraryCounts;
}

#####################
# INPUT SUBROUTINES #
#####################

sub readReadRegions {
    my $readRegionFile = shift;
    my @readRegions;
    open(RRF,$readRegionFile);
    while(<RRF>) {
        chomp;
	unless(/\#/) {	  
	    my($id,$location,$length) = split(/\t/,$_);
	    push(@readRegions,[$id,$location,$length]);
	}
    }
    return \@readRegions;
}

sub loadReadGffList {
    my($readGffListFile,$parameters) = @_;
    my(%dataSet,%sampleTotal,@sampleList,%hitCount);
    my $maxHitCount = $parameters->{maxHitCount} or die "FAIL: could not load maxHitCount\n";
    my $thisSource="";
    my $total=0;
    open(RFILE,$readGffListFile) or die "could not open $readGffListFile\n";
    while(<RFILE>) {
        chomp;
        my($sample,$gffFile)=split;
	my %allHitCount;
        open(GFILE,$gffFile) or die "could not open $gffFile\n";
        while(<GFILE>) {
            unless(/\#/) {
                chomp;
                my($chrom,$source,$side,$start,$stop,$score,$strand,$phase,$info)=split(/\t/,$_);
		my $readId;
		if(($readId) = $info =~ /ID=(.*?);/) {
		    $allHitCount{$readId}++;
		} elsif(($readId) = $info =~ /Match (.*)/) {
		    $allHitCount{$readId}++;
		}
	    }
	}
	seek(GFILE,0,0);
        while(<GFILE>) {
            unless(/\#/) {
                chomp;
                my($chrom,$source,$side,$start,$stop,$score,$strand,$phase,$info)=split(/\t/,$_);
		my $readId;
		if(/ID=(.*?);/) {
		    ($readId) = $info =~ /ID=(.*?);/;
		} elsif(/Match .*/) {
		    ($readId) = $info =~ /Match (.*)/;
		}
		unless($allHitCount{$readId} > $maxHitCount) {
		    # store the sample name after visiting the first valid read.
		    unless($sampleTotal{$sample}) {			
			push(@sampleList,$sample);
		    }
		    if(my($count) = $readId =~ /.*_x(\d+)/) {
			# this could be improved. could include the copy number in dataSet
			for(my $c=0;$c<$count;$c++) {
			    push(@{$dataSet{$chrom}{$strand}},[$start,$stop,$sample,$readId]);
			}
			$sampleTotal{$sample} += $count;
			# each line of this gff is only one hit to the genome, 
			# even though maybe sequenced in multiple copies.
			$hitCount{$sample.$readId}++;
		    } else {
			push(@{$dataSet{$chrom}{$strand}},[$start,$stop,$sample,$readId]);
			$sampleTotal{$sample}++;
			$hitCount{$sample.$readId}++;
		    }
		}
	    }
	}
	close(GFILE);
    }
    close(RFILE);
    return(\%dataSet,\%sampleTotal,\@sampleList,\%hitCount);
}

sub readRNAFoldOutput {
    my $rnaFoldOutput = shift;
    my %folds;
    my($id,$sequence,$tag,$location);
    open(MFILE,$rnaFoldOutput);
    while(<MFILE>) {
        chomp;
        if(/>/) {
            s/>//;
            ($tag,$location) = split;
        } elsif(/[ACGU]/) {
            $sequence = $_;
        } elsif(/(.*)\s+\((-\d+\.\d+)\)/) {
            my $fold = $1;
            my $mfe = $2;
            my($chrom,$start,$stop,$strand) = parseLocation($location);
            push(@{$folds{$chrom}{$strand}},[$start,$stop,$tag,$mfe,$sequence,$fold]);
            if(length($sequence) != length($fold)) {
                die "lengths don't match for $id $tag\n";
            }
        }
    }
    return \%folds;
}

sub readMirRegionsFile {
    my $mirRegionsFile = shift;
    my %folds;
    my($id,$sequence,$tag,$location);
    open(MFILE,$mirRegionsFile);
    while(<MFILE>) {
        chomp;
	my($tag,$location,$leftCenter,$rightCenter,$mfe,$sequence,$fold) = split;
	my($chrom,$start,$stop,$strand) = parseLocation($location);
	push(@{$folds{$chrom}{$strand}},[$tag,$start,$stop,$leftCenter,$rightCenter,$mfe,$sequence,$fold]);
    }
    close(MFILE);
    return \%folds;
}

sub readRepeatRegions {
    my $repeatRegionFile = shift;
    my %repeatRegions;
    open(RRF,$repeatRegionFile) or die "could not open $repeatRegionFile\n";
    while(<RRF>) {
        chomp;
        my($chrom,$start,$stop) = split(/\t/,$_);
        push(@{$repeatRegions{$chrom}},[$start,$stop]);
    }
    return \%repeatRegions;
}

sub readtRNAScanFile {
    my $tRNAScanFile = shift;
    my %tRNAs;
    my $RECORD = 0;
    if($tRNAScanFile) {
        open(TRNA,$tRNAScanFile);
        while(<TRNA>) {
            chomp;
            my($id,$count,$tStart,$tStop,$type,$antiCodon,$seqId,$length,$score) = split(/\s+/);
            if($RECORD) {
                $tRNAs{$id} = 1;
            }
            if($id =~ /\-+/) {
                $RECORD = 1;
            }
        }
    }
    return \%tRNAs;
}

sub readFoldFile {
    my $foldFile = shift;
    my @foldInfo;
    open(FF,$foldFile) or die "could not open $foldFile\n";
    while(<FF>) {
        chomp;
        my($name,$location,$mfe,$seq,$fold) = split(/\t/);
        push(@foldInfo,[$name,$location,$mfe,$seq,$fold]);
    }
    return \@foldInfo;
}

sub readAllDistinctReadsFile {
    my($allDistinctReadsFile) = @_;
    my($chromTitle,$chromStrand,$startTitle,$stopTitle,$totalTitle,$totalHitCountTitle,$readTotalTitle,@libraryNames);
    my(%allDistinctReads,%libraryNames);
    open(ADR,$allDistinctReadsFile);
    while(<ADR>) {
	chomp;
	if(/\#/) {
	    ($chromTitle,$chromStrand,$startTitle,$stopTitle,$totalTitle,$totalHitCountTitle,$readTotalTitle,@libraryNames)=split(/\t/,$_);
	} else {
	    my($chrom,$strand,$dStart,$dStop,$dTotal,$dTotalHitCount,$dReadTotal,@libraryCounts) = split(/\t/);
	    # note changing the order here:
	    my %libraryCounts;
	    for(my $i=0;$i<@libraryNames;$i++) {
		$libraryCounts{$libraryNames[$i]} = $libraryCounts[$i];
	    }
	    push(@{$allDistinctReads{$chrom}{$strand}},[$dStart,$dStop,$dTotal,\%libraryCounts,$dTotalHitCount,$dReadTotal]);
	}
    }
    return(\%allDistinctReads,\@libraryNames);
}

sub readDistinctReadsFile {
    my $mirDistinctReadsFile = shift;
    open(MDR,$mirDistinctReadsFile) or die "could not open $mirDistinctReadsFile\n";
    my($nameTitle,$startPosTitle,$stopPosTitle,@libraryNames);
    my %distinctReadsHash;
    while(<MDR>) {
        chomp;
        if(/\#/) {
           ($nameTitle,$startPosTitle,$stopPosTitle,@libraryNames)=split(/\t/,$_);
       } else {
           my($name,$startPos,$stopPos,@libraryCounts)=split(/\t/,$_);
	   my %libraryCountHash;
	   my $totalCount = 0;
	   for(my $i=0;$i<@libraryCounts;$i++) {
	       $totalCount += $libraryCounts[$i];
	       $libraryCountHash{$libraryNames[$i]} = $libraryCounts[$i];
	   }
           push(@{$distinctReadsHash{$name}},[$startPos,$stopPos,
					      $totalCount,\%libraryCountHash]);
       }
    }
    close(MDR);
    foreach my $name (keys %distinctReadsHash) {
	# sort by abundance highest to lowest...
	@{$distinctReadsHash{$name}} = sort {$b->[2] <=> $a->[2]} @{$distinctReadsHash{$name}}
    }
    return(\%distinctReadsHash,\@libraryNames);
}

sub readProductFile {
    my $productFile = shift;
    open(PF,$productFile) or die "could not open $productFile\n";
    my $COUNT = 1;
    my @productList;
    while(<PF>) {
	chomp;	
	my($tag,$side,$type,$prodCount,$maxProdCount,$maxStart,$maxStop,$sequence,@samples) = split(/\t/);
	my $id = $tag . "_" . $side . "_" . $type;
	$COUNT++;
	push(@productList,[$id,$prodCount,$maxProdCount,$maxStart,$maxStop,$sequence,\@samples]);
    }
    return \@productList;
}

sub readHairpinListFile {
    my $hairpinFile = shift;
    my %hairpinList;
    open(MFILE,$hairpinFile);
    while(<MFILE>) {
        chomp;
        unless(/\#/) {
            my($tag,$chrom,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,
	       $arpd,$trpd,$urf,$ahc,$afh,$pbp,$shifted) = split(/\t/);
	    push(@{$hairpinList{$chrom}{$strand}},[$tag,$start,$stop,$hStart,$hStop,$total,$mfe,
						   $arpd,$trpd,$urf,$ahc,$afh,$pbp,$shifted]);
        }
    }
    return \%hairpinList;
}

sub readMirListFile {
    my $mirListFile = shift;
    my %mirList;
    open(MFILE,$mirListFile);
    while(<MFILE>) {
        chomp;
        unless(/\#/) {
            my($tag,$chrom,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,$sequence,$fold) = split(/\t/);
	    push(@{$mirList{$chrom}{$strand}},[$tag,$start,$stop,$hStart,$hStop,$total,$mfe,$sequence,$fold]);
        }
    }
    return \%mirList;
}

sub readMirProductFile {
    my $mirProductFile = shift;
    my %mirProducts;
    my($tagT,$sideT,$typeT,$prodCountT,$maxProdT,$relStartT,$relStopT,$seqT,@libraries);
    open(HPF,$mirProductFile) or die "could not open $mirProductFile\n";
    while(<HPF>) {
	chomp;
	unless(/\#/) {
	    my($tag,$side,$type,$prodCount,$maxProdCount,$relStart,$relStop,$productSequence,@sampleCounts) = split(/\t/);
	    push(@{$mirProducts{$tag}},[$side,$type,$prodCount,$maxProdCount,$relStart,$relStop,$productSequence,@sampleCounts]);
	} else {
	    ($tagT,$sideT,$typeT,$prodCountT,$maxProdT,$relStartT,$relStopT,$seqT,@libraries) = split(/\t/);
	}
    }
    return(\%mirProducts,\@libraries);
}

1;
__END__

=head1 NAME

miRTRAP (miR Tests for Read Analysis and Prediction) - Perl module (and script) for discovering microRNAs from highthroughput sequencing data.

=head1 SYNOPSIS

#the following code is the contents of the associated "printReadRegions.pl" script.
use miRTRAP;
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

=head1 DESCRIPTION

MicroRNAs (miRs) have been broadly implicated in animal development and disease. However, the systematic, whole-genome identification of miRs is complicated by their small size.  Current identification methods rely mainly on the projected stability of putative stem-loop structures (pre-miRs) and on sequence comparisons with known miRs. 

This method, miRTRAP, incorporates the mechanisms of miR biogenesis, and includes additional criteria regarding the prevalence and quality of small RNAs arising from the antisense strand and neighboring loci.

Please see associate scripts and README file for more information and examples.

=head2 EXPORT

None by default.

=head1 SEE ALSO

http://flybuzz.berkeley.edu/miRTRAP.html

=head1 AUTHOR

David Hendrix, davidhendrix@berkeley.edu

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by David Hendrix

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=cut
