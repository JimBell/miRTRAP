package miRTRAP1_6;
use Bio::DB::Sam;
use Memory::Usage;
use TimeKeeper;
use RNA;
use strict;
use warnings;


sub loadDefaultParameters {
    my $parameters = {
	"aapdLimit" => 3,
	"ahcLimit" => 5,
	"afhLimit" => 0.5,
	"minSameShift" => 7,
	"minBothShift" => 7,
	"NEIGHBOR" => 1,
	"CHECK_EXONS" => 0,
	"nonMirNeighborLimit" => 20,
	"neighborWindow" => 1000,
	"minGap" => 3,
	"minMirCount" => 1,  
	"MaxLength" => 160,
	"lengthMin" => 20,
	"totalLength" => 150,
	"distanceMin" => 10,
	"CountMax" => 5,
	"HitCountMax" => 50,
	"MaxMajor" => 22,
	"minMajor" => -44,
	"reverseMax" => 0.05,
	"countMinLocus" => 1,   
	"fivePrimeHetMax" => 0.5,
	"shiftMin" => 7,
	"OverlapMin" => 2,
	"BPdensityLimit" => 0.6,
	"hairpinShortLength" => 12,
	"InHairpinBuffer" => 3,
	"OutHairpinBuffer" => 3,
	"RangeOfHairpin" => 70,
	"outputPrefix" => "readRegions",
	"bamListFile" => "",
	"RepeatRegionsFile" => "",
	"genomeDir" => "",
	"SizesFile" => "",
	"geneModels" => "",
	"validList" => "",
	"knownMirsFile" => "",
	"MaxShortProductSize" => "14",
	"minLongProductSize" => "28",
	"LoadFromConfigFile" => "",
	"minOverlapAbundanceLimit" => 0,  #used in determinine whether to find a same shift, both shift, or aapd score for overlapping producs
	"maxSplitProductAbundance" => 0.05,
	"maxOverlapProductAbundance" => 0.05
    };
    return $parameters;
}

sub createOutputFileParameters {
    my($parameters) = shift;

    my $filePrefix = $parameters->{outputPrefix} or die "FAIL: filePrefix not loaded.\n";
    $parameters->{distancesFile} = $filePrefix."_distances.txt";
    $parameters->{readRegionsFile} = $filePrefix.".txt";
    $parameters->{longReadRegionsFile} = $filePrefix."_long.txt";
    $parameters->{allReadRegionsFile} = $filePrefix."_full.txt";
    $parameters->{readRegionsFastaFile} = $filePrefix . ".fasta";
    $parameters->{tRNAScanFasta} = $filePrefix . "_trnaScanFasta.fasta";
    $parameters->{tRNAScanOutputFile} = $filePrefix . ".trna";
    $parameters->{filteredCandidatesFile} = $filePrefix . "_filteredCandidates.txt";
    $parameters->{hairpinsFile} = $filePrefix . "_hairpins.txt";
    $parameters->{distinctReadsFile} = $filePrefix . "_distinctReads.txt";
    $parameters->{productFile} = $filePrefix . "_products.txt";    
    return $parameters;
}

sub readConfigFile {
    my($configFile,$parameters) = @_;
    open(CONFIGFILE,$configFile) or die "FAIL: could not open $configFile\n";
    while(<CONFIGFILE>) {
	chomp;
	my($key,$value) = split(/\s+=\s+/);
	$parameters->{$key} = $value;
    }
    return $parameters;
}


######################################
# SEQUENCE / GENOME TOOLS            #
######################################

sub mapBamStrand {
    my $strand = shift;
    if($strand eq "-1") {
	return "-";
    }
    return "+";
}

sub revStrand {
    my $strand = shift;
    if($strand eq "-") {
	return "+";
    }
    return "-";
}


sub reverseComplement {
# Returns the reverse complement of the input sequence.
    my($seq)=@_;
    $seq =~ tr/acgtuACGTU/tgcaaTGCAA/;
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

    my %sequences;
    my $tag;

    open(GFILE,$genomeFile) or die "could not open $genomeFile\n";
    while(<GFILE>) {
	s/\r?\n|\r/\n/g; # in case of dos or mac
        chomp;
        if(/>/) {
            s/>//g;
            my @terms = split;
	    $tag = shift(@terms);
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
	# converting from 1-based to 0-based below:
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

sub readChromLengths {
    my($chromSizesFile) = @_;
    my %chromLengths;
    open(CSF, "$chromSizesFile");
    while (<CSF>) {
	chomp;
	my($chrom,$size) = split(/\s+/);
	$chromLengths{$chrom} = $size;
    }
    close(CSF);
    return \%chromLengths;
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

sub isInt {
    my $num = shift;
    if ($num =~ /^\d+$/) {
	return 1;
    }
    return 0;
}

#########################
# FOLD PROCESSING TOOLS #
#########################

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

sub getHairpinCenters {
    my($fold,$parameters) = @_;
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
    my @centers = getHairpinCenters($fold,$parameters);
    my($mergedCenters,$MERGED) = mergeCenters(\@centers,$fold,$parameters);
    while($MERGED) {
	($mergedCenters,$MERGED) = mergeCenters($mergedCenters,$fold,$parameters);
    }
    return @{$mergedCenters};
}

sub mergeCenters {
    my($centers,$fold,$parameters) = @_;
    my $minLength = $parameters->{lengthMin};
    my $shortLength = $parameters->{hairpinShortLength};
    my @mergedCenters;
    my $MERGED = 0;
    for(my $i=0;$i<@{$centers};$i++) {
	#print "checking $center\n";
	my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$centers->[$i]);
	my($leftCenter, $rightCenter) = @{$centers->[$i]};
	my $leftLength = $leftCenter - $leftEnd + 1;
	my $rightLength = $rightEnd - $rightCenter + 1;
	if(($leftLength >= $minLength)&&($rightLength >= $shortLength)) {
	    # sufficiently large of a loop, keep as is.
	    push(@mergedCenters,$centers->[$i]);
	} elsif(($leftLength >= $minLength)&&($rightLength < $shortLength)) {
	    # potentially merge the loop.
	    if($i<@{$centers}-1) {
		# if there is a neighboring center
		my($leftEnd2,$rightEnd2) = extendHairpinToFullEnds($fold,$centers->[$i+1]);
		my($leftCenter2, $rightCenter2) = @{$centers->[$i+1]};
		my $leftLength2 = $leftCenter2 - $leftEnd2 + 1;
		my $rightLength2 = $rightEnd2 - $rightCenter2 + 1;
		if($leftLength2 < $shortLength) {
		    #print "merging $center $leftCenter2,$rightCenter2\n";
		    my($tightLeftEnd,$tightRightEnd) = extendHairpinToEnds($fold,$centers->[$i]);
		    my($tightLeftEnd2,$tightRightEnd2) = extendHairpinToEnds($fold,$centers->[$i+1]);
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
		my($leftEnd2,$rightEnd2) = extendHairpinToFullEnds($fold,$centers->[$i-1]);
		my($leftCenter2, $rightCenter2) = @{$centers->[$i-1]};
		my $leftLength2 = $leftCenter2 - $leftEnd2 + 1;
		my $rightLength2 = $rightEnd2 - $rightCenter2 + 1;
		if($rightLength2 < $shortLength) {
		    #print "merging $center $leftCenter2,$rightCenter2\n";
		    my($tightLeftEnd,$tightRightEnd) = extendHairpinToEnds($fold,$centers->[$i]);
		    my($tightLeftEnd2,$tightRightEnd2) = extendHairpinToEnds($fold,$centers->[$i-1]);
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

sub getFullHairpinWindow {
    # this method returns the window of paired bases, extends beyond minor hairpins    
    my($fold,$center) = @_;;
    return extendHairpinToFullEnds($fold,$center);
}

sub getBasePairs {
    my($center,$fold) = @_;
    my($leftCenter, $rightCenter) = @{$center};
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
        die "crapped out at $center:\n$fold\n";
    }
}

sub getMaxHairpinLength {
    my($fold,$center) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$center);
    my $leftLength = $leftCenter - $leftEnd + 1;
    my $rightLength = $rightEnd - $rightCenter + 1;
    return max($leftLength,$rightLength);
}

sub extendHairpinToEnds {
    my($fold,$center) = @_;
    my($leftCenter, $rightCenter) = @{$center};
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
    my($fold,$center) = @_;
    my($leftCenter, $rightCenter) = @{$center};
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
        die "extendHairpinFullToEnds FAIL: ($center)\n$fold\n$L\n$R\n\n";	
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


#############################
# BAM RETRIEVAL SUBROUTINES #
#############################

sub loadBamFile {
    my $bamFile = shift;
    my $bam = Bio::DB::Bam->open( $bamFile );
    return $bam;
}

sub loadBamIndex {
    my $bamFile = shift;
    my $reIndex;  #changed to 1 if the index file doesn't exist
    my $bamIndex =  Bio::DB::Bam->index($bamFile,$reIndex);
    die "failed to load index for $bamFile\n" if ($reIndex);
    return $bamIndex;
}

sub loadBamList {
    my $bamListFile = shift;
    open(BLF,$bamListFile) or die "could not open $bamListFile\n";
    my @bamList;
    while(<BLF>) {
	chomp;
	my($label,$bamFile) = split;
	my $bam = loadBamFile($bamFile);
	my $bamHeader = $bam->header;
	my $bamIndex = loadBamIndex($bamFile);
	my $totalMapped = getBamTotal($bamFile);
	push(@bamList,[$bam,$bamHeader,$bamIndex,$bamFile,$label,$totalMapped]);
    }
    return \@bamList;
}

sub getBamTotal {
    my $bamFile = shift;
    my $bamFlagstatFile = $bamFile . ".flagstat";
    unless(-e $bamFlagstatFile) {
	# if flagstat file hasn't been created, make one.
	system("samtools flagstat $bamFile > $bamFlagstatFile");
    }    
    open(BFS,$bamFlagstatFile) or die "Could not open $bamFlagstatFile.\n";
    while(<BFS>) {
	if(/(\d+) \+ \d+ mapped/) {
	    close(BFS);
	    return $1;
	}
    }
    close(BFS);
    # if we made it here, something went wrong with parsing the flagstat file
    die "Error. Could not parse the flagstat file here: $bamFlagstatFile. Older version of samtools perhaps?\n"
}

sub retrieveReadData {
    my($location,$bamList) = @_;
    my($chrom,$start,$stop,$strand) = miRTRAP1_6::parseLocation($location);
    $location = "$chrom:$start-$stop";  
    my %distinctReadCount;  #normalized by dividing by the hitCount for each read before totalling the number of reads
    my %libraryCounts;  #normalized by dividing by the hitCount for each read before totalling the number of reads
    my %adjustedLibraryCounts;
    my %hitCounts;
    my %distinctReads;
    my %readTotal;
    my %adjustedSeqCount;
    my @sampleList;
    my %uniqueReadCount;
    my %readCount;
    my %adjustedReadCount;   #normalized by dividing by the hitCount
    foreach my $bamData (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamData};
	my($tid,$chromStart,$chromStop) = $bamHeader->parse_region("$chrom:$start-$stop"); #converts start and end to zero based coordinates
	push(@sampleList, $sample);
	my $callBack = sub {
	    my($alignment,$data) = @_;
	    my($chromStart,$chromStop,$distinctReadCount,$libraryCounts,$adjustedLibraryCounts,$hitCounts,
	       $adjustedSeqCount,$readTotal,$uniqueReadCount,$readCount,$adjustedReadCount) = @{$data};
	    my $id = $alignment->qname;
	    my $rStart = $alignment->start;  #returned in 1 based coordinates
	    my $rStop = $alignment->end;  #returned in 1 based coordinates
	    my $seq = $alignment->qseq;
	    my $hitCount = $alignment->get_tag_values('NH');
	    my $rStrand = miRTRAP1_6::mapBamStrand($alignment->strand);
	    my $count;
	    if ($id =~ /.*_x(\d+)$/) {
		($count) = $id =~ /.*_x(\d+)$/;
	    } else {
		$count = 1;
	    }
	    $seq = miRTRAP1_6::reverseComplement($seq) if ($strand eq '-');
	    my $length = length($seq);
	    if(($chromStart <= $rStart)&&($rStop <= $chromStop)) {
		$distinctReadCount->{$rStrand}{$rStart}{$seq} += $count;
		$libraryCounts->{$rStrand}{$rStart}{$seq}{$sample} += $count;
		$adjustedLibraryCounts->{$rStrand}{$rStart}{$seq}{$sample} += $count / $hitCount;
		$hitCounts->{$seq} += $count * $hitCount;
		$adjustedSeqCount->{$seq} += $count/$hitCount;
		$readTotal->{$seq} += $count;
		$uniqueReadCount->{$rStrand} += $count if ($hitCount == 1);
		$readCount->{$rStrand} += $count;
		$adjustedReadCount->{$rStrand} += $count / $hitCount;
	    }
	};
	# end callBack
	my $callBackData = [$chromStart,$chromStop,\%distinctReadCount,\%libraryCounts,\%adjustedLibraryCounts,\%hitCounts,
			    \%adjustedSeqCount,\%readTotal,\%uniqueReadCount,\%readCount,\%adjustedReadCount];
	my $code = $bamIndex->fetch($bam,$tid,$chromStart,$chromStop,$callBack,$callBackData);
    }
    foreach my $rStrand (keys %distinctReadCount) {
	$uniqueReadCount{$rStrand} = 0 unless ($uniqueReadCount{$rStrand});
	$readCount{$rStrand} = 0 unless ($readCount{$rStrand});
	$adjustedReadCount{$rStrand} = 0 unless ($adjustedReadCount{$rStrand});
	foreach my $rStart (keys %{$distinctReadCount{$rStrand}}) {
	    foreach my $seq (keys %{$distinctReadCount{$rStrand}{$rStart}}) {
		foreach my $sample (@sampleList) {
		    unless ($libraryCounts{$rStrand}{$rStart}{$seq}{$sample}) {
			$libraryCounts{$rStrand}{$rStart}{$seq}{$sample} = 0;
			$adjustedLibraryCounts{$rStrand}{$rStart}{$seq}{$sample} = 0;
		    }
		}
		my $avgHitCount = $hitCounts{$seq} / $readTotal{$seq};
		my $adjustedSeqCount = $adjustedSeqCount{$seq} / $readTotal{$seq};
		push(@{$distinctReads{$rStrand}}, [$rStart, $rStart + length($seq) - 1, $distinctReadCount{$rStrand}{$rStart}{$seq}, $avgHitCount, 
						   $libraryCounts{$rStrand}{$rStart}{$seq}, $seq, $adjustedSeqCount, 
						   $adjustedLibraryCounts{$rStrand}{$rStart}{$seq}]);
	    }
	}
    }
    # sort distinct reads by abundance.  This is important in getProductInfo
    foreach my $strand (keys %distinctReads) {
	@{$distinctReads{$strand}} = sort {$b->[2] <=> $a->[2]} @{$distinctReads{$strand}};
    }
    return \%distinctReads, \%uniqueReadCount, \%readCount, \%adjustedReadCount;
}


sub getTagValue {
    my $tags = shift;
    my $tagName = shift;
    my $tagType = shift;
    foreach my $tag (@{$tags}) {
	my($tagValue) = $tag =~ /^$tagName:$tagType:(.*)/;
	if ($tagValue) {
	    return $tagValue;
	}
    }
    return 0;
}

sub getDistinctReadCounts {
    my $location = shift;
    my $bamList = shift;
    my %readCounts;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    foreach my $bamEntry (@{$bamList}) {
	my($label, $bam, $total, $bamFile) = @{$bamEntry};	
	my @reads = $bam->get_features_by_location(-seq_id => $chrom, -start  => $start, -end => $stop);
	foreach my $read (@reads) {
	    my $rStart = $read->start;
	    my $rStop = $read->end;
	    # 5' end methdod.
	    my $rStrand = mapBamStrand($read->strand);
	    # we could remove this restriction later for requiring the same strand
	    if($strand eq $rStrand) {
		$readCounts{$rStart."x".$rStop}{$label}++;
	    }
	}
    }
    my @readCounts;
    foreach my $key (keys %readCounts) {
	my($rStart,$rStop) = $key =~ /(\d+)x(\d+)/;
	my @counts;
	my $totalCounts = 0;
	foreach my $bamEntry (@{$bamList}) {
	    my($label, $bam, $total, $bamFile) = @{$bamEntry};
	    my $count = $readCounts{$key}{$label} ? $readCounts{$key}{$label} : 0;
	    $totalCounts += $count;
	    push(@counts,$count);
	}
	push(@readCounts,[$rStart,$rStop,$totalCounts,\@counts]);
    }
    return \@readCounts;
}

###############################
# tRNA PROCESSING SUBROUTINES #
###############################

sub printTRNAFastaFile {
    my($hairpinsFile,$tRNAScanFastaFile) = @_;
    open(HP,"$hairpinsFile") or die "failed to open $hairpinsFile for reading in printTRNAFastaFile()";
    open(TRNAFF, ">$tRNAScanFastaFile") or die "failed to open $tRNAScanFastaFile for writing in printTRNAFastaFile()";
    while(<HP>) {
	my($id,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$adjTotalReads,$mfe,$seq,$fold) = split(/\t/);
	my $location = "$chrom:$start..$stop:$strand";
	print TRNAFF ">$id\t$location\n$seq\n";
    }
    close(HP);
    close(TRNAFF);
}  


###############################
# READ PROCESSING SUBROUTINES #
###############################

sub printReadRegions {
    my $bamList = shift;
    my $chromLengths = shift;
    my $repeatRegions = shift;
    my $parameters = shift;
    my $maxLength = $parameters->{MaxLength} or die "FAIL: maxLength not loaded (not found in parameters.)\n";
    my $readRegionsFile = $parameters->{readRegionsFile} or die "FAIL: readRegionsFile not loaded (not found in parameters.)\n";
    my $longReadRegionsFile = $parameters->{longReadRegionsFile} or die "FAIL: longReadRegionsFile not loaded (not found in parameters.)\n";
    my $allReadRegionsFile = $parameters->{allReadRegionsFile} or die "FAIL: allReadRegionsFile not loaded (not found in parameters.)(\n";
    my $initWindowLength = 200000;
    my $shortReadsCount = 1;     # for regions shorter than 100bp.
    my $longReadsCount = 1;     # for regions longer than 100bp.
    my $repeatReadsCount = 1;     # for repeat regions.

    my $mu = Memory::Usage->new();

    open(RF,">".$readRegionsFile) or die "Failed to load $readRegionsFile for writing\n";
    open(LRF,">".$longReadRegionsFile) or die "Failed to load $longReadRegionsFile for writing\n";
    open(ARF,">".$allReadRegionsFile) or die "Failed to load $allReadRegionsFile for writing\n";
    foreach my $chrom (keys %{$chromLengths}) {
	my $start = 0;
	while($start <= $chromLengths->{$chrom}) {
	    # push out the windowLength until no reads extend beyond the ends.
	    my $windowLength = $initWindowLength <= $chromLengths->{$chrom} ? $initWindowLength : $chromLengths->{$chrom};
	    # now the window is defined so that all contiguous read regions are properly contained within it.
	    my($plusArray,$minusArray) = getReadCountArrays($bamList,$chrom,$start,$start+$windowLength);
	    # we now have plus and minus arrays for thie window.
 	    my $posReadRegions = extractReadRegionList($plusArray,$start,$windowLength,"+");
 	    my $negReadRegions = extractReadRegionList($minusArray,$start,$windowLength,"-");
	    # sort readRegions in descending order by rrStart position
	    my @readRegions = sort {$b->[0] <=> $a->[0]} (@{$posReadRegions}, @{$negReadRegions});
	    my $newStart = $start+$windowLength;

	    $mu->record("readRegions obtained for $chrom:$start-$newStart");

	    foreach my $readRegion (@readRegions) {
		my($rrStart,$rrStop,$rrStrand,$rrMaxCount) = @{$readRegion};
		my $prrStart = $rrStart-1; # printable 0-based start for bed file
		if($rrStop < $newStart) {
		    unless(insideList($rrStart,$rrStop,$repeatRegions->{$chrom})) {
			# print the info.
			my $length = $rrStop - $rrStart;
			if($length <= $maxLength) {
			    print ARF "$chrom\t$prrStart\t$rrStop\trr$shortReadsCount\t$rrMaxCount\t$rrStrand\n";
			    print RF "$chrom\t$prrStart\t$rrStop\trr$shortReadsCount\t$rrMaxCount\t$rrStrand\n";
			    $shortReadsCount++;
			} else {
			    print ARF "$chrom\t$prrStart\t$rrStop\tlrr$longReadsCount\t$rrMaxCount\t$rrStrand\tRejected:Region too long\n";
			    print LRF "$chrom\t$prrStart\t$rrStop\tlrr$longReadsCount\t$rrMaxCount\t$rrStrand\n";
			    $longReadsCount++;
			}
		    } else {
			print ARF "$chrom\t$prrStart\t$rrStop\trrr$repeatReadsCount\t$rrMaxCount\t$rrStrand\tRejected: region within repeat.\n";
			$repeatReadsCount++;
		    }
		} else {
		    $newStart = $rrStart;
		}
	    }
	    # increment the start to beyond this window.
	    die "Huge read region at $chrom:$start.." . ($start+$windowLength) . "\nPlease report to miRTRAP developers.\n" if($start == $newStart);
	    # put the next window at beginning of last read region
	    $start = $newStart;
	}
    }
    close(RF);
    close(LRF);
    close(ARF);
    $mu->dump();
}

sub getReadCountArrays {
    my($bamList,$chrom,$windowStart,$windowStop) = @_;
    my $location = "$chrom:$windowStart-$windowStop";
    my(@plusArray,@minusArray);
    foreach my $bamElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
	my @bamARGS = ("samtools", "view", $bamFile, $location);
	open(BAM, "-|", @bamARGS) or die "Failed to open pipe to @bamARGS";
	while (<BAM>) {
	    chomp;
	    unless( /^@/ ) {
		my($id, $flag, $refSeqName, $rStart, $mapq, $cigar, $mrnm, $mpos, $iSize, $seq, $qualityScore, @tags) = split(/\t/);
		my $readCount;
		if ($id =~ /.*_x(\d+)$/) {
		    ($readCount) = $id =~ /.*_x(\d+)$/;
		} else {
		    $readCount = 1;
		}
		my $rStrand = $flag & 16 ? '-': '+';  #bitwise anding to see if sequence is reverse complemented
		my $length = length($seq);
		my $rStop = $rStart + $length - 1;
		if($rStart >= $windowStart) {
		    if($rStrand eq "+") {
			# in plus strand
			for(my $i=$rStart;$i<=$rStop;$i++) {
			    $plusArray[$i-$windowStart]++;     #generating arrays relative to windowStart
			}
		    } else {
			# in minus strand
			for(my $i=$rStart;$i<=$rStop;$i++) {
			    $minusArray[$i-$windowStart]++;
			}
		    }
		} else {
		    die "Found read at $chrom:$windowStart..$windowStop that extends beyond beginning: $chrom:$rStart..$rStop:$rStrand\n";
		}
	    }
	}
    }
    return(\@plusArray,\@minusArray);
}

sub extractReadRegionList {
    my($readCountArray,$start,$windowLength,$strand) = @_;
    my $regionStart;
    my $regionStop;
    my @readRegions;
    my $INREGION = 0;
    my $maxCount = 0;
    for (my $j=0; $j<@{$readCountArray}; $j++) {
	if($readCountArray->[$j]) {
	    if($INREGION) {
		# extend the regionStop
		$regionStop = $j;
	    } else {
		# assume the region is only 1bp long initially.
		$regionStart = $j;
		$regionStop = $j;
		$INREGION = 1;
	    }
	    $maxCount = $readCountArray->[$j] if($maxCount < $readCountArray->[$j]);
	} else {
	    # no reads here, store regionious region if just leaving one.
	    if($INREGION) {
		push(@readRegions, [$regionStart+$start,$regionStop+$start,$strand,$maxCount]);
		$INREGION = 0;
		$maxCount = 0;
	    }
	}
    }
    # if still in region at the end of the loop push start and stop into readRegions array
    if($INREGION) {
	push(@readRegions, [$regionStart+$start,$regionStop+$start,$strand,$maxCount]);
	$INREGION = 0;
    }
    
    return \@readRegions;
}

sub extendSequenceToFillBuffer {
    my($chrom, $start, $stop, $strand, $chromLength, $parameters) = @_;
    my $bufferLength = getBufferLength($start,$stop,$parameters);
    my $newStart = $start-$bufferLength > 1 ? $start - $bufferLength : 1;
    my $newStop = $stop+$bufferLength <= $chromLength ? $stop+$bufferLength : $chromLength;
    my $newLocation = $chrom . ":" . $newStart . ".." . $newStop  . ":" . $strand;
    return $newLocation;
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
    my($bamList,$readRegionsFile,$genomeDir,$chromLengths,$parameters) = @_;
    my $testTime = 0;
    my $testMemory = 0;
    my $testReadCounts = 0;
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    my $readsLessThanCountMinLocus = 0; #number of hairpins with fewer than count min locus reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noCenterCount = 0;  #number of hairpins without centers (hopefully this will stay 0)
    my $mu = Memory::Usage->new() if ($testMemory);
    my $tk = TimeKeeper->new() if ($testTime);
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    open(FRC,">".$filteredCandidatesFile) or die "failed to open $filteredCandidatesFile for writing\n";
    close(FRC);
    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
    open(HPL,">".$productFile) or die "failed to open $productFile for writing";
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tafh\tvalid\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";

    $tk->start("TotalTime") if ($testTime);

    open(RRF, $readRegionsFile) or die "failed to open $readRegionsFile\n";
    my $prevChrom;
    my $genome;
    while (<RRF>) {
        chomp;
	unless(/\#/) {	  
	    my($chrom,$start,$stop,$id,$score,$strand) = split(/\t/,$_);
	    $start++; #converting to one based
	    #the chroms in RRF should be sorted at this point so the next line of code will only be executed once per chrom
	    if ($prevChrom) {
		$genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	    } else {
		$genome = loadGenome("$genomeDir/$chrom.fa");	
	    }
	    $prevChrom = $chrom;
	    my $location = extendSequenceToFillBuffer($chrom,$start,$stop,$strand,$chromLengths->{$chrom},$parameters);
	    ($chrom,$start,$stop,$strand) = parseLocation($location);
	    my $sequence = getSequence($location,$genome);
    	    $tk->start("retrieveReadData") if ($testTime);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList);
	    $tk->end("retrieveReadData") if ($testTime);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		if ($adjustedReadCount->{$strand} < $parameters->{countMinLocus}) {
		    $readsLessThanCountMinLocus++;
		    print "$id has fewer reads than count min locus\n" if ($testReadCounts);
		}
		$readsCount++;
		$tk->start("buildProducts") if ($testTime);
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		$tk->end("buildProducts") if ($testTime);
		$tk->start("getFold") if ($testTime);
		my($fold,$mfe) = getFold($chrom,$sequence,$senseProducts,$id);
		$tk->end("getFold") if ($testTime);
		$tk->start("getMergedHairpinCenters") if ($testTime);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		$tk->end("getMergedHairpinCenters") if ($testTime);
		$centersReadsCount += scalar(@centers);
		unless (scalar(@centers)) {
		    $noCenterCount++;
		    print "$id has no centers\n" if ($testReadCounts);
		    print "$fold\n" if ($testReadCounts);
		}
		my $COUNT= 0;
		foreach my $center (@centers) {
		    my $asciiVal = ord('a') + $COUNT;
		    my $label = chr($asciiVal);
		    my $newId = $id . $label;
		    $COUNT++;
		    $tk->start("getBasePairs") if ($testTime);
		    my $basePairs = getBasePairs($center,$fold);
		    $tk->end("getBasePairs") if ($testTime);
		    $tk->start("getMaxHairpinLength") if ($testTime);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    $tk->end("getMaxHairpinLength") if ($testTime);
		    $tk->start("getFullHairpinWindow") if ($testTime);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    $tk->end("getFullHairpinWindow") if ($testTime);
		    if($hairpinLength >= $minLength) {
			$tk->start("getProductInfo") if ($testTime);
			my $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
			$tk->end("getProductInfo") if ($testTime);
			$tk->start("rebuildProducts") if ($testTime);
			$senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
			$tk->end("rebuildProducts") if ($testTime);
			$senseProducts=addProductSize($senseProducts,$parameters);
			if(goodProducts($senseProducts,$parameters)) {
			    my $revStrand = revStrand($strand);
			    $tk->start("buildProducts") if ($testTime);
			    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
			    $tk->end("buildProducts") if ($testTime);
			    $tk->start("getProductInfo") if ($testTime);
			    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
			    $tk->end("getProductInfo") if ($testTime);
			    $tk->start("rebuildProducts") if ($testTime);
			    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
			    $tk->end("rebuildProducts") if ($testTime);
			    $antiSenseProducts=addProductSize($antiSenseProducts,$parameters);
			    $tk->start("getAdjTotalProductReads") if ($testTime);
			    my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
			    $tk->end("getAdjTotalProductReads") if ($testTime);
			    $tk->start("getAdjTotalProductReads") if ($testTime);
			    my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
			    $tk->end("getAdjTotalProductReads") if ($testTime);
			    my($PASS,$REASON) = plausibleReads($adjTotalProductReads,$adjTotalRevProductReads,$senseProducts,$antiSenseProducts,
							       $parameters);
			    if($PASS) {
				$tk->start("getReverseProductDisplacement") if ($testTime);
				my($tpd,$totalRP)=newGetReverseProductDisplacement($senseProducts,
										$antiSenseProducts,
										$adjTotalProductReads,
										$parameters,$newId);
				$tk->end("getReverseProductDisplacement") if ($testTime);
				my $apd = $totalRP ? $tpd/$totalRP : 0.0;
				my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
				$tk->start("computeMaxProdHitCount") if ($testTime);
				my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
				$tk->end("computeMaxProdHitCount") if ($testTime);
				$tk->start("computeMaxProdFivePrimeHet") if ($testTime);
				my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
				$tk->end("computeMaxProdFivePrimeHet") if ($testTime);
				$tk->start("computeProductBasePairing") if ($testTime);
				my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
				$tk->end("computeProductBasePairing") if ($testTime);
				$tk->start("computeMaxSameShift") if ($testTime);
				my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
				$tk->end("computeMaxSameShift") if ($testTime);
				$tk->start("computeMaxBothShift") if ($testTime);
				my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters,$newId);
				$tk->end("computeMaxBothShift") if ($testTime);
				my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
				my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
				# total is the read count for the whole hairpin
				my $totalSense = 0;
				foreach my $product (@{$senseProducts}) {
				    my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
				       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				    my $length = $relStop - $relStart + 1;
				    my $productSequence = substr($sequence,$relStart,$length);
				    $totalSense += $adjProdCount;
				    $tk->start("getAdjTotalLibraryCounts") if ($testTime);
				    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
				    $tk->end("getAdjTotalLibraryCounts") if ($testTime);
				    print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
				    print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
				    print HPL "$productSequence";
				    foreach my $bamElement (@{$bamList}) {
					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
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
					    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					    printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
					}
					print HDR "\n";
				    }				
				}
				#############################
				# ANTISENSE PRODUCTS
				my $totalAntisense = 0;
				foreach my $product (@{$antiSenseProducts}) {
				    my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
				       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				    $newType = "as-".$newType;
				    $totalAntisense += $adjProdCount;
				    my $length = $relStop - $relStart + 1;
				    my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
				    $tk->start("getAdjTotalLibraryCounts") if ($testTime);
				    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
				    $tk->end("getAdjTotalLibraryCounts") if ($testTime);
				    print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
				    print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
				    print HPL "$productSequence";
				    foreach my $bamElement (@{$bamList}) {
					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
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
					    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
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
				printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\n",
				       $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$splitProductAbundance,$overlapProductAbundance);
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
		    } 
		}
	    } else {
		$noReadsCount++;
		print "no reads found at $id\t$location\n" if ($testReadCounts);
	    }
	}
    }
    $mu->record("after finishing") if ($testMemory);
    $tk->end("TotalTime") if ($testTime);
    $tk->dumpTimes() if ($testTime);

    print "hairpins with reads = $readsCount\n" if ($testReadCounts);
    print "hairpins without reads = $noReadsCount\n" if ($testReadCounts);
    print "hairpin centers = $centersReadsCount\n" if ($testReadCounts);
    print "hairpins without centers (not including those with no reads) = $noCenterCount\n" if ($testReadCounts);
    print "hairpins with reads less than countMinLocus (not including those with no reads) = $readsLessThanCountMinLocus\n" if ($testReadCounts);

    close(RRF);
    close(HPL);
    close(HDR);
    $mu->dump() if ($testMemory);
}


sub mirPreprocess {
    my($bamList,$hairpins,$genomeDir,$chromLengths,$parameters) = @_;
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    open(FRC,">".$filteredCandidatesFile) or die "failed to open $filteredCandidatesFile for writing\n";
    close(FRC);
    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
    open(HPL,">".$productFile) or die "failed to open $productFile for writing";
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tafh\tvalid\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";

    my $prevChrom;
    my $genome;
    foreach my $chrom (keys %{$hairpins}) {
	if ($prevChrom) {
	    $genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	} else {
	    $genome = loadGenome("$genomeDir/$chrom.fa");	
	}
	$prevChrom = $chrom;
	foreach my $hairpin (@{$hairpins->{$chrom}}) {
	    my ($start,$stop,$strand,$id,$mirName) = @$hairpin;
	    my $location = "$chrom:$start..$stop:$strand";
	    my $sequence = getSequence($location,$genome);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		my($fold,$mfe) = getFold($chrom,$sequence,$senseProducts,$mirName);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		my $COUNT= 0;
		foreach my $center (@centers) {
#		    my $asciiVal = ord('a') + $COUNT;
#		    my $label = chr($asciiVal);
#		    my $newId = $mirName . $label;
		    $COUNT++;
		    my $basePairs = getBasePairs($center,$fold);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    my $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
		    $senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
		    $senseProducts=addProductSize($senseProducts,$parameters);
		    my $revStrand = revStrand($strand);
		    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
		    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
		    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
		    $antiSenseProducts=addProductSize($antiSenseProducts,$parameters);
		    my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
		    my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
		    my($tpd,$totalRP)=newGetReverseProductDisplacement($senseProducts,
										$antiSenseProducts,
										$adjTotalProductReads,
										$parameters,$mirName);
		    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
		    my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
		    my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
		    my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
		    my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
		    my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
		    my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters,$mirName);
		    my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
		    my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
		    # total is the read count for the whole hairpin
		    my $totalSense = 0;
		    foreach my $product (@{$senseProducts}) {
			my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			my $length = $relStop - $relStart + 1;
			my $productSequence = substr($sequence,$relStart,$length);
			$totalSense += $adjProdCount;
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			print HPL "$mirName\t$side\t$newType\t$adjProdCount\t";
			print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$prodList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
			    print HDR "$offset\t$relStart\t$relStop\t$seq";
			    foreach my $bamElement (@{$bamList}) {
				my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
			    }
			    print HDR "\n";
			}				
		    }
		    #############################
		    # ANTISENSE PRODUCTS
		    my $totalAntisense = 0;
		    foreach my $product (@{$antiSenseProducts}) {
			my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			$newType = "as-".$newType;
			$totalAntisense += $adjProdCount;
			my $length = $relStop - $relStart + 1;
			my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			print HPL "$mirName\t$side\t$newType\t$adjProdCount\t";
			print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$prodList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
			    print HDR "$offset\t$relStart\t$relStop\t$seq";
			    foreach my $bamElement (@{$bamList}) {
				my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
			    }
			    print HDR "\n";
			}				
		    }
		    #############################
		    #############################
		    my($leftCenter, $rightCenter) = @{$center};
		    print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
		    print HL "$sequence\t$fold\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\n",
			   $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$splitProductAbundance,$overlapProductAbundance);

		}
	    } else {
		print "no reads found at $mirName\t$location\n";
	    }
	}
    }

    close(RRF);
    close(HPL);
    close(HDR);
}



sub oldmirPreprocess {
    my($bamList,$hairpins,$genomeDir,$chromLengths,$parameters) = @_;
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    open(FRC,">".$filteredCandidatesFile) or die "failed to open $filteredCandidatesFile for writing\n";
    close(FRC);
    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
    open(HPL,">".$productFile) or die "failed to open $productFile for writing";
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tafh\tvalid\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";

    my $prevChrom;
    my $genome;
    foreach my $chrom (keys %{$hairpins}) {
	if ($prevChrom) {
	    $genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	} else {
	    $genome = loadGenome("$genomeDir/$chrom.fa");	
	}
	$prevChrom = $chrom;
	foreach my $hairpin (@{$hairpins->{$chrom}}) {
	    my ($start, $stop, $strand, $id, $mirName) = @$hairpin;
	    my $location = "$chrom:$start..$stop:$strand";
	    my $sequence = getSequence($location,$genome);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		$readsCount++;
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		my($fold,$mfe) = getFold($chrom,$sequence,$senseProducts,$mirName);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		$centersReadsCount += scalar(@centers);
		my $COUNT= 0;
		foreach my $center (@centers) {
		    my $asciiVal = ord('a') + $COUNT;
		    my $label = chr($asciiVal);
		    my $newId = $mirName . $label;
		    print ("$mirName has more than one label\n") if ($label eq 'b');
		    $COUNT++;
		    my $basePairs = getBasePairs($center,$fold);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    my $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
		    $senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
		    my $revStrand = revStrand($strand);
		    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
		    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
		    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
		    my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
		    my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
		    my($tpd,$totalRP)=getReverseProductDisplacement($senseProducts,
								    $antiSenseProducts,
								    $adjTotalProductReads,
								    $parameters);
		    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
		    my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
		    my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
		    my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
		    my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
		    my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
		    my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters);
		    my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
		    my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
		    # total is the read count for the whole hairpin
		    my $totalSense = 0;
		    foreach my $product (@{$senseProducts}) {
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
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
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
				    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				    printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
					}
				print HDR "\n";
			    }				
		    }
		    #############################
		    # ANTISENSE PRODUCTS
		    my $totalAntisense = 0;
		    foreach my $product (@{$antiSenseProducts}) {
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
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
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
				my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
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
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\n",
			   $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$splitProductAbundance,$overlapProductAbundance);
		}
	    } else {
		$noReadsCount++;
		print "no reads found at $mirName\t$location\n";
	    }
	}
    }
    print "hairpins with reads = $readsCount\n";
    print "hairpin centers with reads = $centersReadsCount\n";
    print "hairpins without reads = $noReadsCount\n";

    close(RRF);
    close(HPL);
    close(HDR);
}


sub getFold {
    my($chrom,$sequence,$senseProducts,$hairpinId) = @_;
    my $test = 0;
    my $majorProduct = $senseProducts->[0];  #products must be sorted by abundance at this point
    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
    print "$hairpinId\n" if ($test);
    my $productSequence = substr($sequence,$maxRelStart,$maxRelStop-$maxRelStart+1);
    print "Most Abundant Product Sequence:$productSequence\n" if ($test);
    my $productSeqLength = length($productSequence);
    my $maskedSequence = maskSequence($sequence,$maxRelStart,$maxRelStop);
    print "sequence:\t\t$sequence\n" if ($test);
    print "masked sequence:\t$maskedSequence\n" if ($test);
    my ($productDuplexCoords,$seqDuplexCoords) = getDuplex($productSequence,$maskedSequence);
    print "Most Abundant Product Coords: $maxRelStart $maxRelStop\n" if ($test);
    print "productDuplexCoords: $productDuplexCoords->[0] $productDuplexCoords->[1]\n" if ($test);
    print "sequenceDuplexCoords: $seqDuplexCoords->[0] $seqDuplexCoords->[1]\n" if ($test);
    print "duplexSequence: " . substr($maskedSequence,$seqDuplexCoords->[0],$seqDuplexCoords->[1]-$seqDuplexCoords->[0]+1) . "\n" if ($test);

    my $duplexSide = ($maxRelStart > $seqDuplexCoords->[0]) ? '5p' : '3p'; #in reality this is only true for real mirs but it works for our purpose here
    print "duplex side = $duplexSide\n" if ($test);

    
    #mirStarShift is the length of the mir* not covered by the duplex
    #in rare cases the duplex region may overlap a mor but at the moment there is no way to check if a product is a mor until the hairpin is folded
    my $mirStarShift = getMirStarShift($seqDuplexCoords,$duplexSide,$senseProducts);
    print "mirStarShift = $mirStarShift\n" if ($test);
    #mirShift is the buffer to add due to the duplexFold not finding the duplex energy for the full product sequence
    my $mirShift = ($duplexSide eq '5p') ? $productSeqLength - $productDuplexCoords->[1] - 1 : $productDuplexCoords->[0];
    print "mirShift = $mirShift\n" if ($test);
    my $maxMirShift = ($mirStarShift > $mirShift) ? $mirStarShift: $mirShift;
    my($foldStart,$foldStop) = extendHairpinToCoverMirs($majorProduct,$seqDuplexCoords,$maxMirShift,$duplexSide);
    print "adding $maxMirShift to both sides and 2 to the $duplexSide side\nfoldStart = $foldStart foldStop = $foldStop\n" if ($test);
    $foldStart = 0 if ($foldStart < 0);
    $foldStop = length($sequence) - 1 if ($foldStop > (length($sequence)-1));

    my $foldSequence = substr($sequence,$foldStart,$foldStop-$foldStart+1);
    print "foldSequence=\t$foldSequence\n" if ($test);
    my($fold,$mfe) = RNA::fold($foldSequence);
    my $leftOutRegion = '+' x $foldStart;
    my $rightOutRegion = '+' x (length($sequence)-$foldStop-1);
    my $newFold = $leftOutRegion . $fold . $rightOutRegion;
    print "sequence=\t$sequence\n" if ($test);
    print "fold=    \t$newFold\n" if ($test);
    print "\n\n" if ($test);
    return $newFold,$mfe,$seqDuplexCoords;
}

sub getDuplex {
    my($seq1,$seq2) = @_;
    my $dup = RNA::duplexfold($seq1,$seq2);
    my $dupStructure = $dup->{structure};
    my($seq1DupStruct,$seq2DupStruct) = $dupStructure =~ /^(.*)&(.*)/;  #needed to get the lengths of the duplex coords.
    my $seq1Length = length($seq1);
    my $seq1DupStructLength = length($seq1DupStruct);             
    my $seq2DupStructLength = length($seq2DupStruct);                          
    my @seq1DuplexCoords = ($dup->{i}-$seq1DupStructLength, $dup->{i} - 1);   
    my @seq2DuplexCoords = ($dup->{j} - 1,$dup->{j}+$seq2DupStructLength-2); 
    return (\@seq1DuplexCoords, \@seq2DuplexCoords);
}

sub getMirStarShift {
    my($duplexCoords,$duplexSide,$products) = @_;
    my($leftDuplexCoord,$rightDuplexCoord) = @{$duplexCoords};
    if ($duplexSide eq '5p') {
	foreach my $product (@{$products}) {
	    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$product};
	    if (getOverlap($maxRelStart,$maxRelStop,$leftDuplexCoord,$rightDuplexCoord)) {
		if ($maxRelStart < $leftDuplexCoord) {
		    return $leftDuplexCoord - $maxRelStart;
		}
	    }
	}
    } elsif ($duplexSide eq '3p') {
	foreach my $product (@{$products}) {
	    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$product};
	    if (getOverlap($maxRelStart,$maxRelStop,$leftDuplexCoord,$rightDuplexCoord)) {
		if ($maxRelStop > $rightDuplexCoord) {
		    return $maxRelStop - $rightDuplexCoord;
		}
	    }
	}
    } else {
	die "no duplexSide given\n";
    }
    return 0;
}

sub extendHairpinToCoverMirs {
    my($majorProduct,$seqDuplexCoords,$mirBuffer,$mirStarSide) = @_;
    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
    my($foldStart,$foldStop);
    if ($mirStarSide eq '5p') {
	$foldStart = $seqDuplexCoords->[0] - $mirBuffer;
	$foldStop = $maxRelStop + $mirBuffer;
	$foldStart = $foldStart - 2;  #since there are two unpaired nucleotides on the hairpin
    } else {
	$foldStart = $maxRelStart - $mirBuffer;
	$foldStop = $seqDuplexCoords->[1] + $mirBuffer;
	$foldStop = $foldStop + 2;   #since there are two unpaired nucleotides on the hairpin
    }

    return ($foldStart,$foldStop);
}

sub oldGetFold {
    my($chrom,$sequence,$senseProducts,$hairpinId) = @_;
    my $test = 1;
    my $majorProduct = $senseProducts->[0];  #products must be sorted by this point
    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
    print "$hairpinId\n" if ($test);
    print "maxRelStart = $maxRelStart\tmaxRelStop = $maxRelStop\n" if ($test);
    my $productSequence = substr($sequence,$maxRelStart,$maxRelStop-$maxRelStart+1); #converted to zero based to get substring
    print "product Sequence:$productSequence\n" if ($test);
    my $maskedSequence = maskSequence($sequence,$maxRelStart,$maxRelStop);
    print "sequence:\t\t$sequence\n" if ($test);
    print "masked sequence:\t$maskedSequence\n" if ($test);
    my $dup = RNA::duplexfold($productSequence,$maskedSequence);
    my $dupStructure = $dup->{structure};
    my($productSeqDupStruct,$seqDupStruct) = $dupStructure =~ /^(.*)&(.*)/;  #needed to get the lengths.
    my $productSeqLength = $maxRelStop - $maxRelStart + 1;
    my $productSeqDupLength = length($productSeqDupStruct);             
    my $seqDupLength = length($seqDupStruct);                           #a better way is to use the lengths from the structure to find coords
    my @productDuplexCoords = ($dup->{i}-$productSeqDupLength, $dup->{i} - 1);   #duplex coords of main product in 0 based
    print "ProductCoords: $maxRelStart $maxRelStop\n" if ($test);
    print "productDuplexCoords: $productDuplexCoords[0] $productDuplexCoords[1]\n" if ($test);
    my @seqDuplexCoords = ($dup->{j} - 1,$dup->{j}+$seqDupLength-2); #duplex coords of sequence in 0 based
    print "sequenceDuplexCoords: $seqDuplexCoords[0] $seqDuplexCoords[1]\n" if ($test);
    print "duplexSequence: " . substr($maskedSequence,$seqDuplexCoords[0],$seqDuplexCoords[1]-$seqDuplexCoords[0]+1) . "\n" if ($test);
    my $pairedDuplexSide = ($maxRelStart > $seqDuplexCoords[0]) ? 'left' : 'right';
    print "duplex side = $pairedDuplexSide\n" if ($test);
    #overlappingProductBuffer is the buffer to add to the region that is due to overlapping products
    my $overlappingProductBuffer = getOverlappingProductBuffer(\@seqDuplexCoords,$pairedDuplexSide,$senseProducts);
    print "overlappingProductBuffer = $overlappingProductBuffer\n" if ($test);
    #maxProductBuffer is the buffer to add due to the duplexFold not finding the duplex energy for the full product sequence
    my $maxProductBuffer = ($pairedDuplexSide eq 'left') ? $productSeqLength - $productDuplexCoords[1] - 1 : $productDuplexCoords[0];
    print "maxProductBuffer = $maxProductBuffer\n" if ($test);

    my($foldStart,$foldStop,$foldSequence);    
    if ($pairedDuplexSide eq 'left') {
	if ($overlappingProductBuffer) {
	    $foldStart = $seqDuplexCoords[0] - $overlappingProductBuffer - 2;
	    print "subtracting overlapping product buffer: FoldStart = $seqDuplexCoords[0] - $overlappingProductBuffer - 2 = $foldStart\n" if ($test);
	} else {
	    $foldStart = $seqDuplexCoords[0] - $maxProductBuffer - 2;
	    print "no Overlapping product: subtracting maxProductBuffer: FoldStart = $seqDuplexCoords[0] - $maxProductBuffer - 2 = $foldStart\n" if ($test);
	}
	$foldStop = $maxRelStop + $maxProductBuffer;
	print "adding maxProductBuffer: FoldStop = $maxRelStop + $maxProductBuffer = $foldStop\n" if ($test);
    } elsif ($pairedDuplexSide eq 'right') {
	if ($overlappingProductBuffer) {
	    $foldStop = $seqDuplexCoords[1] + $overlappingProductBuffer + 2;
	    print "adding overlapping product buffer: FoldStop = $seqDuplexCoords[1] + $overlappingProductBuffer + 2 = $foldStop\n" if ($test);
	} else {
	    $foldStop = $seqDuplexCoords[1] + $maxProductBuffer + 2;
	    print "no Overlapping product: adding maxProductBuffer: FoldStop = $seqDuplexCoords[1] + $maxProductBuffer + 2 = $foldStop\n" if ($test);
	}
	$foldStart = $maxRelStart - $maxProductBuffer;
	print "subtracting maxProductBuffer: FoldStart = $maxRelStart - $maxProductBuffer = $foldStart\n" if ($test);
    }
    $foldStart = 0 if ($foldStart < 0);
    $foldStop = length($sequence) - 1 if ($foldStop > (length($sequence)-1));
    $foldSequence = substr($sequence,$foldStart,$foldStop-$foldStart+1);
    print "foldSequence=\t$foldSequence\n" if ($test);
    my($fold,$mfe) = RNA::fold($foldSequence);
    my $leftOutRegion = '+' x $foldStart;
    my $rightOutRegion = '+' x (length($sequence)-$foldStop-1);
    my $newFold = $leftOutRegion . $fold . $rightOutRegion;
    print "sequence=\t$sequence\n" if ($test);
    print "fold=    \t$newFold\n" if ($test);
    print "\n\n" if ($test);
    return $newFold,$mfe,\@seqDuplexCoords;
}


sub computeOverlapProductAbundance {
    my($productInfo,$parameters) = @_;
    my $minOverlap = $parameters->{OverlapMin};
    my $productTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProdCount;
    }
    my @sortedProducts = sort {$a->[3] <=> $b->[3]} @{$productInfo};  #sorted from least to most abundant
    my $maxOverlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1) = @{$product1};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @sortedProducts; $j++) {
		my $product2 = $sortedProducts[$j];
		my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2) = @{$product2};
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			if ($overlap > $minOverlap) {
			    my $overlapProductAbundance = $adjProdCount1 / $productTotal;
			    $maxOverlapProductAbundance = $overlapProductAbundance if ($overlapProductAbundance > $maxOverlapProductAbundance);
			}
		    }
		}
	    }
	}	
    }
    return $maxOverlapProductAbundance;
}

#this version of computeOverlapProductAbundance does not take into account product size
sub oldComputeOverlapProductAbundance {
    my($productInfo,$parameters) = @_;
    my $productTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProdCount;
    }
    my @sortedProducts = sort {$a->[3] <=> $b->[3]} @{$productInfo};  #sorted from least to most abundant
    my $maxOverlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1) = @{$product1};
	for (my $j = $i + 1; $j < @sortedProducts; $j++) {
	    my $product2 = $sortedProducts[$j];
	    my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2) = @{$product2};
	    if (getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
		my $overlapProductAbundance = $adjProdCount1 / $productTotal;
		$maxOverlapProductAbundance = $overlapProductAbundance if ($overlapProductAbundance > $maxOverlapProductAbundance);
	    }
	}	
    }
    return $maxOverlapProductAbundance;
}

sub computeSplitProductAbundance {
    my($productInfo,$parameters) = @_;
    my $productTotal = 0;
    my $splitProductTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	unless(($newType =~ /long/) || ($newType =~ /short/)) {
	    $splitProductTotal += $adjProdCount if ($newType eq 'split');
	    $productTotal += $adjProdCount;
	}
    }
    return ($productTotal == 0) ? 0 : $splitProductTotal / $productTotal;
}

#this version of computeSplitProductAbundance does not take into account product size
sub oldComputeSplitProductAbundance {
    my($productInfo,$parameters) = @_;
    my $productTotal = 0;
    my $splitProductTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$splitProductTotal += $adjProdCount if ($newType eq 'split');
	$productTotal += $adjProdCount;
    }
    return ($productTotal == 0) ? 0 : $splitProductTotal / $productTotal;
}

#sub oldGetFold {
#    my($chrom,$sequence,$senseProducts,$hairpinId) = @_;
#    my $majorProduct = $senseProducts->[0];  #products must be sorted by this point
#    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
#    print "$hairpinId\t$productReads\n";
#    my $productSequence = substr($sequence,$maxRelStart-1,$maxRelStop-$maxRelStart+1); #converted to zero based to get substring
#    print "product Sequence:\n$productSequence\n";
#    my $maskedSequence = maskSequence($sequence,$maxRelStart,$maxRelStop);
#    print "masked sequence:\n$maskedSequence\n";
#    my $dup = RNA::duplexfold($productSequence,$maskedSequence);
#    my $dupStructure = $dup->{structure};
#    my($productSeqDupStruct,$seqDupStruct) = $dupStructure =~ /^(.*)&(.*)/;  #needed to get the lengths.
#    my $productSeqLength = length($productSeqDupStruct);              #example in RNA.pm for duplexFold was incorrect in how to determine the duplex coords.
#    my $seqLength = length($seqDupStruct);                           #a better way is to use the lengths from the structure to find coords
#    my @productDuplexCoords = ($dup->{i}-$productSeqLength, $dup->{i}-1);  #these coordinants are in 0 based
#    print "productDuplexCoords: $productDuplexCoords[0] $productDuplexCoords[1]\n";
#    my @seqDuplexCoords = ($dup->{j}-1,$dup->{j}+$seqLength-2);
#    print "sequenceDuplexCoords: $seqDuplexCoords[0] $seqDuplexCoords[1]\n";
#    my($foldStart,$foldStop,$foldSequence);
#    my $buffer = $productDuplexCoords[0];  #buffer is added if the duplex found isn't exact
#    print "buffer: $buffer\n";
    #in the following if statement maxRelStart and maxRelStop are one based.  foldStart and foldStop will be zero based.
#    if ($maxRelStart < $seqDuplexCoords[0]) {
#	print "maxRelStart = $maxRelStart < $seqDuplexCoords[0] = seqDuplexCoords[0]\n";
	#the part of the sequence the product duplexes with is after it
#	$foldStart = ($maxRelStart <= $buffer + 1) ? 0 : $maxRelStart - $buffer - 1;  
#	$foldStop = ($seqDuplexCoords[1] >= length($sequence) - $buffer - 1) ? $seqDuplexCoords[1] : $seqDuplexCoords[1] + $buffer;
#	print "foldStart = $foldStart\tfoldStop=$foldStop\n";
#   } else {
#       print "maxRelStart = $maxRelStart > $seqDuplexCoords[0] = seqDuplexCoords[0]\n";
       #the part of the sequence the product duplexes with is before it
#       $foldStart = ($seqDuplexCoords[0] <= $buffer) ? 0 : $seqDuplexCoords[0] - $buffer;
#       $foldStop = ($maxRelStop >= length($sequence) - $buffer) ? $maxRelStop - 1 : $maxRelStop + $buffer - 1;
#	print "foldStart = $foldStart\tfoldStop=$foldStop\n";
#    }
#    $foldSequence = substr($sequence,$foldStart,$foldStop-$foldStart+1);
#    print "foldSequence = $foldSequence\n";
#    my($fold,$mfe) = RNA::fold($foldSequence);
#    my $leftOutRegion = '+' x $foldStart;
#    my $rightOutRegion = '+' x (length($sequence)-$foldStop-1);
#    my $newFold = $leftOutRegion . $fold . $rightOutRegion;
#    return $newFold,$mfe,\@seqDuplexCoords;
#}


sub buildProducts {
    my($location,$distinctReads,$productStrand,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $PRODUCTCOUNT = 1;
    my %productHash;
    my %relStartHash;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts) = @{$dRead};
	my $relStart = ($strand eq "+") ? $dStart - $start : $stop - $dStop;
	my $relStop = ($strand eq "+") ? $dStop - $start : $stop - $dStart;
	#offsets will be determined later after the fold
	my $parsedRead = [$relStart,$relStop,$dStart,$total,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts,$productStrand];
	if (my $id = overlapCurrent(\%productHash,$parsedRead,$parameters)) {
	    push(@{$productHash{$id}},$parsedRead);
	} else {
	    push(@{$productHash{$PRODUCTCOUNT}},$parsedRead);
	    $PRODUCTCOUNT++;
	}
    }
    my @productInfo;
    foreach my $id (keys %productHash) {
	my $adjTotal=0;
	my($maxRelStart,$maxRelStop,$maxGStart);
	my $FIRST=1;
	my $adjProductCount = 0;
	my $adjMaxProductCount = 0;
	my %dReadCounts;
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$storedRead};
	    $dReadCounts{$relStart."x".$relStop} += $readCount * $adjustedSeqCount;
	    $relStartHash{$relStart."x".$relStop} = $relStart;
	}
	my @keys = sort {$dReadCounts{$b} <=> $dReadCounts{$a} || $relStartHash{$a} <=> $relStartHash{$b}} keys %dReadCounts;
#	my @keys = sort {$dReadCounts{$b} <=> $dReadCounts{$a}} keys %dReadCounts;
	my $maxKey = shift(@keys);
	($maxRelStart,$maxRelStop) = split(/x/,$maxKey); 
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$storedRead};
	    $adjProductCount += $readCount * $adjustedSeqCount;
	    if(($relStart == $maxRelStart)&&($relStop == $maxRelStop)) {
		$adjMaxProductCount += $readCount * $adjustedSeqCount;
		$maxGStart = $dStart;
	    }
	}
	push(@productInfo,[$id,$productHash{$id},$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand]);
    }
    my @sortedProducts = sort {$b->[3] <=> $a->[3]} @productInfo;
    return \@sortedProducts;
}

sub maskSequence {
    my($sequence,$start,$stop) = @_;
    my $seqLength = length($sequence);
    my $maskLength = $stop-$start+1;
    my $maskedSequence;
    #note the product sequence is in one base but the sequence is in 0 base
    if ($stop == $seqLength) {
	return substr($sequence,0,$start) . 'N' x ($maskLength - 1);
    }
    return substr($sequence,0,$start) . 'N' x ($maskLength) . substr($sequence, $stop + 1, $seqLength - $stop - 1);
}

sub extractProducts {
    my($center,$fold,$basePairs,$location,$distinctReads,$drStrand,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $productInfo = getProductInfo($center,$fold,$basePairs,
				     $start,$stop,$strand,
				     $distinctReads,$drStrand,$parameters);
    $productInfo = rebuildProducts($center,$fold,$basePairs,
				    $productInfo,$parameters);    
    return $productInfo;
}

sub computeMaxSameShift {
    my($location,$strandDistinctReads,$productInfo, $parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $newParameters = {};
    foreach my $key (keys %{$parameters}) {
	$newParameters->{$key} = $parameters->{$key};
    }
    $newParameters->{distanceMin} = $parameters->{shiftMin};
    my $adjTotalReads = getAdjTotalReads($strandDistinctReads,$parameters);
    return getSameArmShift($productInfo,$adjTotalReads,$newParameters);
}

sub computeMaxBothShift {
    my($basePairs,$location,$strandDistinctReads,$productInfo,$parameters,$mirId) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $newParameters = {};
    foreach my $key (keys %{$parameters}) {
	$newParameters->{$key} = $parameters->{$key};
    }
    $newParameters->{distanceMin} = $parameters->{shiftMin};
    my $adjTotalReads = getAdjTotalReads($strandDistinctReads,$parameters);
    if(bothArmProducts($basePairs,$productInfo,$newParameters)) {
	#print "possible both arm shift\n";
	return newGetBothArmShift($basePairs,$productInfo,$adjTotalReads,$newParameters,$location,$mirId);
    }
    return 0;
}

sub bothArmProducts {
    my($basePairs,$productInfo,$parameters) = @_;
    my %sideHash;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless(($newType =~ /long/) || ($newType =~ /short/)) {
	     if($side eq "3p") {
		 if(overlapsBasePairs($basePairs,$relStart,$relStop)) {
		     $sideHash{$side}++;
		 }
	     } else {
		 $sideHash{$side}++;
	     }
	}
    }
    if(($sideHash{"5p"})&&($sideHash{"3p"})) {
	return 1;
    }
    return 0;
}

sub oldBothArmProducts {
    my($basePairs,$productInfo,$parameters) = @_;
    my %sideHash;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
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


sub newGetBothArmShift {
    my($basePairs,$productInfo,$adjTotalReads,$parameters,$location,$newId) = @_;
    my $test = 0;
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my %FIVE_PRIME_USED;
    my %THREE_PRIME_USED;
    my @productOverlaps;
    my $maxShift = 0;
    print "\n\n$newId\n" if ($test);
    print "minReadFraction = $minReadFraction\tcountThreshold = $countThreshold\tminOverlap = $minOverlap\n" if ($test);
     for(my $i=0;$i<@{$productInfo};$i++) {
	 my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	     if($adjProdCount1 >= $countThreshold) {	
		 if($side1 eq "5p") {
		     print "adjProdCount1 = $adjProdCount1 >= $countThreshold = countThreshold for $side1-$newType1\n" if ($test); 
		     for(my $j=0;$j<@{$productInfo};$j++) {
			 if($i != $j) {
			     my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			     unless(($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				 if($adjProdCount2 >= $countThreshold) {
				     if($side2 eq "3p") {				
					 print "\tadjProdCount2 = $adjProdCount2 >= $countThreshold = countThreshold for $side2-$newType2\n" if ($test);
					 if(overlapsBasePairs($basePairs,$relStart2,$relStop2)) {
					     my($relStart2bp,$relStop2bp) = getOuterBasePairsFivePrime($basePairs,
												       $relStart2,
												       $relStop2,
												       $location);
					     my $overlap = getOverlap($relStart1,$relStop1,
								      $relStart2bp,$relStop2bp);
					     my $shift = getShift($relStart1,$relStop1,
								  $relStart2bp,$relStop2bp);
					     push(@productOverlaps, [$overlap,$shift,$i,$j]);
					     print "\tBase pairs for $side1-$newType1 and $side2-$newType2 overlap by $overlap" if ($test);
					     print "(shift = $shift)\n" if ($test);
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
    #finding the maximum shift for products that are sorted from most abundant to least abundant overlap
    @productOverlaps = sort {$b->[0] <=> $a->[0]} @productOverlaps;
    foreach my $item (@productOverlaps) {
	my($overlap,$shift,$i,$j) = @{$item};
	unless(($FIVE_PRIME_USED{$i})||($THREE_PRIME_USED{$j})) {
	    $maxShift = $shift if ($shift > $maxShift);
	    $FIVE_PRIME_USED{$i}++;
	     $THREE_PRIME_USED{$j}++;  
	}
    }
    print "maximum shift =  $maxShift\n" if ($test);
    
    return $maxShift;
}



sub getBothArmShift {
    my($basePairs,$productInfo,$adjTotalReads,$parameters,$location,$newId) = @_;
    my $test = 0;
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $maxShift = 0;
    print "\n\n$newId\n" if ($test);
    print "minReadFraction = $minReadFraction\tcountThreshold = $countThreshold\tminOverlap = $minOverlap\tmaxShfit = $maxShift\n" if ($test);
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    if($adjProdCount1 >= $countThreshold) {	
		if($side1 eq "5p") {
		    print "adjProdCount1 = $adjProdCount1 >= $countThreshold = countThreshold for $side1-$newType1\n" if ($test); 
		    for(my $j=0;$j<@{$productInfo};$j++) {
			if($i != $j) {
			    my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			    unless(($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				if($adjProdCount2 >= $countThreshold) {
				    if($side2 eq "3p") {				
					print "\tadjProdCount2 = $adjProdCount2 >= $countThreshold = countThreshold for $side2-$newType2\n" if ($test);
					if(overlapsBasePairs($basePairs,$relStart2,$relStop2)) {
					    my($relStart2bp,$relStop2bp) = getOuterBasePairsFivePrime($basePairs,
												      $relStart2,
												      $relStop2,
												      $location);
					    my $overlap = getOverlap($relStart1,$relStop1,
								     $relStart2bp,$relStop2bp);
					     print "\tBase pairs for $side1-$newType1 and $side2-$newType2 overlap by $overlap\n" if ($test);
					     my $shift = getShift($relStart1,$relStop1,
								  $relStart2bp,$relStop2bp);
					    print "\tshift = $shift\n" if ($test);
					    if($overlap > $minOverlap) {
						if(abs($shift) > abs($maxShift)) {				
						    $maxShift = $shift;
						    print "\tmaxShift = $maxShift\n" if ($test);
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
	}
    }
    
    return $maxShift;
}

#this version of getBothArmShift does not take into account product size
sub oldGetBothArmShift {
    my($basePairs,$productInfo,$adjTotalReads,$parameters,$location) = @_;
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $maxShift = 0;
    for(my $i=0;$i<@{$productInfo};$i++) {
	 my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,
	    $relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 if($adjProdCount1 >= $countThreshold) {	
	     if($side1 eq "5p") {
		 for(my $j=0;$j<@{$productInfo};$j++) {
		     if($i != $j) {
			 my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,
			    $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			 if($adjProdCount2 >= $countThreshold) {
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
    my($productInfo,$adjTotalReads,$parameters) = @_;
    my $test = 0;
     my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $maxShift = 0;
    print "minOverlap = $minOverlap\tminReadFraction = $minReadFraction\ncountThreshold = $countThreshold\tmaxShift = $maxShift\n" if ($test);
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 unless (($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	     #print "$i: $side1 $newType1 $relStart1 $relStop1 $prodCount1 vs $countThreshold\n";
	     if($adjProdCount1 >= $countThreshold) {
		 print "relStart1 = $relStart1\trelStop1 = $relStop1\n" if ($test);
		 print "adjProdCount1 = $adjProdCount1 >= $countThreshold = countthreshold\n" if ($test);
		 for(my $j=0;$j<@{$productInfo};$j++) {
		     if($i != $j) {
			 my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			 unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
			     if($adjProdCount2 >= $countThreshold) {
				 print "relStart2 = $relStart2\trelStop2 = $relStop2\n" if ($test);
				 print "adjProdCount2 = $adjProdCount2 >= $countThreshold = countThreshold\n" if ($test);
				 my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2);
				 my $shift = getShift($relStart1,$relStop1,$relStart2,$relStop2);
				 print "overlap = $overlap\tshift = $shift\n" if ($test);
				 if($overlap > $minOverlap) {
				     print "overlap = $overlap > $minOverlap = minOverlap\n" if ($test);
				     if($shift > $maxShift) {
					 print "shift = $shift > $maxShift = maxShift\n" if ($test);
					 $maxShift = $shift;
					 print "now maxShift = $shift\n" if ($test);
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }
    }
    print "maxShift = $maxShift\n" if ($test);
    return $maxShift;
}

#this version of getSameArmShift does not take into account product size
sub oldGetSameArmShift {
    my($productInfo,$adjTotalReads,$parameters) = @_;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $maxShift = 0;
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	#print "$i: $side1 $newType1 $relStart1 $relStop1 $prodCount1 vs $countThreshold\n";
	if($adjProdCount1 >= $countThreshold) {	    
	    for(my $j=0;$j<@{$productInfo};$j++) {
		if($i != $j) {
		    my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,
		       $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
		     if($adjProdCount2 >= $countThreshold) {
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
    my($readCount,$revReadCount,$parameters) = @_;
    my $maxReverse = $parameters->{reverseMax};
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
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	my $start = $relStart;
	my $stop = $relStop;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
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


sub newGetReverseProductDisplacement {
    my($productInfo,$revProductInfo,$adjTotalReads,$parameters,$hairpinId) = @_;
    my $test = 0;
    print "\n\n$hairpinId\n" if ($test);
    my $minDispInit = 100;
    # make global variable VVV
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    my %AS_USED;
    my @distances;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    my $start = $relStart;
	    my $stop = $relStop;
	    for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
		my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
		   $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
		unless (($aNewType =~ /long/) || ($aNewType =~ /short/)) {
		    if($aAdjProdCount >= $minFrac*$adjTotalReads) {
			#print "aAdjProdCount = $aAdjProdCount >= $minFrac * $adjTotalReads = minFrac x adjTotalReads\n" if ($test);
			my $aStart = $aRelStart;
			my $aStop = $aRelStop;
			# if the products overlap..
			if((($start <= $aStart)&&($aStart <= $stop))||
			   (($start <= $aStop)&&($aStop <= $stop))||
			   (($aStart <= $start)&&($start <= $aStop))||
			   (($aStart <= $stop)&&($stop <= $aStop))) {
			    my $disp = abs($aStart - $start);
			    push(@distances,[$disp,$i1,$i2]);
			    print "\t$side-$newType overlaps as-$aNewType. Displacement = $disp\n" if ($test);
			}		    
		    }
		}
	    }
	}
    }
    print "\t" if ($test);
    @distances = sort {$a->[0] <=> $b->[0]} @distances;
    foreach my $item (@distances) {
	my($disp,$i1,$i2) = @{$item};
	unless(($USED{$i1})||($AS_USED{$i2})) {
	    $sum += $disp;
	    print "$disp + " if ($test);
	    $total++;
	    $USED{$i1}++;
	    $AS_USED{$i2}++;  
	}
    }
    print " = $sum\n" if ($test);
    if ($total) {
	print "\taapd = $sum / $total = " . $sum / $total . "\n" if ($test);
    }
    return($sum,$total);
}



sub getReverseProductDisplacement {
    my($productInfo,$revProductInfo,$adjTotalReads) = @_;
    my $test = 0;
    my $minDispInit = 100;
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    print "side=$side\ttype=$newType\tadjProdCount=$adjProdCount\tadjMaxProdCount=$adjMaxProdCount\n" if ($test);
	    print "relStart=$relStart\trelStop=$relStop\toffset=$offset\tgStart=$gStart\n" if ($test);
	    my $start = $relStart;
	    my $stop = $relStop;
	    my $minDisp = $minDispInit;
	    print "minDisp=$minDisp\n" if ($test);
	    my $minI2;
	    for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
		my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
		   $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
		unless (($aNewType =~ /long/) || ($aNewType =~ /short/)) {
		    print "aSide=$aType\taType=$aNewType\taAdjProdCount=$aAdjProdCount\taAdjMaxProdCount=$aAdjMaxProdCount\n" if ($test);
		    print "aRelStart=$aRelStart\taRelStop=$aRelStop\taOffset=$aOffset\taGStart=$aGStart\n" if ($test);
		    if($aAdjProdCount >= $minFrac*$adjTotalReads) {
			print "aAdjProdCount = $aAdjProdCount >= $minFrac * $adjTotalReads = minFrac x adjTotalReads\n" if ($test);
			my $aStart = $aRelStart;
			my $aStop = $aRelStop;
			# if the products overlap..
			if((($start <= $aStart)&&($aStart <= $stop))||
			   (($start <= $aStop)&&($aStop <= $stop))||
			   (($aStart <= $start)&&($start <= $aStop))||
			   (($aStart <= $stop)&&($stop <= $aStop))) {
			    my $disp = abs($aStart - $start);
			    print "products overlap: disp = $disp\n" if ($test);
			    unless($USED{$i2}) {
				if($disp < $minDisp) {
				    $minDisp = $disp;
				    $minI2 = $i2;
				    print "disp = $disp < $minDisp = minDisp: now minDisp = $disp\n" if ($test);
				    print "minI2 = $i2\n" if ($test);
				}
			    }
			}		    
		    }
		}
	    }
	    if($minDisp < $minDispInit) {
		$sum += $minDisp;
		$total++;
		$USED{$minI2}++;	    
		print "minDisp = $minDisp < $minDispInit = minDispInit: now sum = sum + $minDisp = $sum\ntotal=$total\n" if ($test);
		print "minI2 = $minI2\n" if ($test);
	    }
	}
    }
    return($sum,$total);
}

#this version of getReverseProductDisplacement does not take into account product size
sub oldGetReverseProductDisplacement {
    my($productInfo,$revProductInfo,$adjTotalReads) = @_;
    my $test = 1;
    my $minDispInit = 100;
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	print "side=$side\ttype=$newType\tadjProdCount=$adjProdCount\tadjMaxProdCount=$adjMaxProdCount\n" if ($test);
	print "relStart=$relStart\trelStop=$relStop\toffset=$offset\tgStart=$gStart\n" if ($test);
	my $start = $relStart;
	my $stop = $relStop;
	my $minDisp = $minDispInit;
	print "minDisp=$minDisp\n" if ($test);
	my $minI2;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
	       $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
	    print "aSide=$aType\taType=$aNewType\taAdjProdCount=$aAdjProdCount\taAdjMaxProdCount=$aAdjMaxProdCount\n" if ($test);
	    print "aRelStart=$aRelStart\taRelStop=$aRelStop\taOffset=$aOffset\taGStart=$aGStart\n" if ($test);
	    if($aAdjProdCount >= $minFrac*$adjTotalReads) {
		print "aAdjProdCount = $aAdjProdCount >= $minFrac * $adjTotalReads = minFrac x adjTotalReads\n" if ($test);
		my $aStart = $aRelStart;
		my $aStop = $aRelStop;
		# if the products overlap..
		if((($start <= $aStart)&&($aStart <= $stop))||
		   (($start <= $aStop)&&($aStop <= $stop))||
		   (($aStart <= $start)&&($start <= $aStop))||
		   (($aStart <= $stop)&&($stop <= $aStop))) {
		    my $disp = abs($aStart - $start);
		    print "products overlap: disp = $disp\n" if ($test);
		    unless($USED{$i2}) {
			if($disp < $minDisp) {
			    $minDisp = $disp;
			    $minI2 = $i2;
			    print "disp = $disp < $minDisp = minDisp: now minDisp = $disp\n" if ($test);
			    print "minI2 = $i2\n" if ($test);
			}
		    }
		}
	    }
	}
	if($minDisp < $minDispInit) {
	    $sum += $minDisp;
	    $total++;
	    $USED{$minI2}++;	    
	    print "minDisp = $minDisp < $minDispInit = minDispInit: now sum = sum + $minDisp = $sum\ntotal=$total\n" if ($test);
	    print "minI2 = $minI2\n" if ($test);
	}
    }
    return($sum,$total);
}

sub plausibleReads {
    my($adjTotalProductReads,$adjTotalRevProductReads,$productInfo,$revProductInfo,$parameters) = @_;
    #    $tk->start("plausibleReads");
    if(doubleStrandedProducts($adjTotalProductReads, $adjTotalRevProductReads, $parameters)) {
	if(overlappingProducts($productInfo,$revProductInfo,$parameters)) {
	    #	    $tk->end("plausibleReads");
	    return (1,"");
	} else {
 #    $tk->end("plausibleReads");
	    return (0,"rejected: has abundant reads on each strand, but none overlap.");
	}
    } else {
	#$tk->end("plausibleReads");
	return (1,"");
    }
}

sub computeProductBasePairing {
    my($center,$productInfo,$basePairs,$parameters) = @_;
    my $minBasePairDensity = 1.0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    if($newType eq "miR") {	    
		my $basePairDensity = computeBasePairDensity($basePairs,$relStart,$relStop);
		if($basePairDensity < $minBasePairDensity) {
		    $minBasePairDensity = $basePairDensity;
		}
	    }
	}
    }
    return $minBasePairDensity;    
}

#this version of computeProductBasePairing does not take into account product size
sub oldComputeProductBasePairing {
    my($center,$productInfo,$basePairs,$parameters) = @_;
    my $minBasePairDensity = 1.0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
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
    #    $tk->start("goodProducts");
    my $maxFivePrimeHet = $parameters->{fivePrimeHetMax} or die "FAIL: no maxFivePrimeHet loaded.\n";
    my $minLocusCount = $parameters->{countMinLocus} or die "FAIL: no minLocusCount loaded.\n";
    my $adjReadCount = getAdjTotalProductReads($productInfo);
    #print "readCount = $readCount\n";
    if($adjReadCount < $minLocusCount) {
	#	$tk->end("goodProducts");
	return 0;
    }
     foreach my $product (@{$productInfo}) {
	 my($side,$newType,$prodList,$prodCount,$maxProdCount) = @{$product};
	 my $fivePrimeHet = computeFivePrimeHet($prodList);
	 #print "prodCount: $prodCount, 5' het: $fivePrimeHet\n";
	 if(($fivePrimeHet < $maxFivePrimeHet)&&($prodCount > 1)) {
	     #	    $tk->end("goodProducts");
	     return 1;
	 }
     }
 #    $tk->end("goodProducts");
    return 0;
}

sub computeMaxProdFivePrimeHet {
    my($productInfo,$parameters) = @_;
    my $maxCount = 0;
    my $maxProdFivePrimeHet = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    if($newType eq "miR") {
		my $fivePrimeHet = computeFivePrimeHet($prodList);
		if($adjProdCount > $maxCount) {
		    $maxProdFivePrimeHet = $fivePrimeHet;
		    $maxCount = $adjMaxProdCount;
		}
	    }
	}
    }
    return $maxProdFivePrimeHet;
}

#this version of computeMaxProdFivePrimeHet does not take into account product size
sub oldComputeMaxProdFivePrimeHet {
    my($productInfo,$parameters) = @_;
    my $maxCount = 0;
    my $maxProdFivePrimeHet = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount) = @{$product};
	if($newType eq "miR") {
	    my $fivePrimeHet = computeFivePrimeHet($prodList);
	    if($adjProdCount > $maxCount) {
		$maxProdFivePrimeHet = $fivePrimeHet;
		$maxCount = $adjMaxProdCount;
	    }
	}
    }
    return $maxProdFivePrimeHet;
}

sub computeFivePrimeHet {
    my($productList) = @_;
    my $FIRST=1;
    my $fivePrimeMaxPos;
    my $adjFivePrimeMaxCount = 0;
    my $adjFivePrimeTotalCount = 0;
    my %startCount;
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts) = @{$read};
	$startCount{$relStart} += $count * $adjustedSeqCount;
    }
    my @starts = sort {$startCount{$b} <=> $startCount{$a}} keys(%startCount);
    my $topStart = shift(@starts);
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts) = @{$read};
	if($relStart == $topStart) {
	    $adjFivePrimeMaxCount += $count * $adjustedSeqCount;
	}
	$adjFivePrimeTotalCount += $count * $adjustedSeqCount;
    }
    my $fivePrimeHet = ($adjFivePrimeTotalCount-$adjFivePrimeMaxCount)/$adjFivePrimeTotalCount;
    return $fivePrimeHet;
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
    my($center,$basePairs,$side,$relStart,$relStop,$newProductInfo,$parameters) = @_;
    my $minShift = $parameters->{shiftMin} or die "FAIL: no minShift loaded.\n";
    foreach my $oProduct (@{$newProductInfo}) {
	my($oSide,$oNewType,$oProductlist,$oAdjustedProductCount,$oAdjustedMaxProdCount,$oRelStart,$oRelStop) = @{$oProduct};
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
    my $inBuffer = $parameters->{InHairpinBuffer};
    my $outBuffer = $parameters->{OutHairpinBuffer};   
    if(($leftEnd-$outBuffer<=$relStart)&&
       ($relStop<=$rightEnd+$outBuffer)) {	    
	return 1;
    }
    return 0;
}

sub onHairpinArm {
    my($center,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my $inBuffer = $parameters->{InHairpinBuffer};
    my $outBuffer = $parameters->{OutHairpinBuffer};   
    if((($leftEnd-$outBuffer<=$relStart)&&($relStop<=$leftCenter+$inBuffer))||
       (($rightCenter-$inBuffer<=$relStart)&&($relStop<=$rightEnd+$outBuffer))) {	    
	return 1;
    }
    return 0;
}

sub closeToHairpin {
    my($center,$relStart,$relStop,$parameters) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my $hairpinRange = $parameters->{RangeOfHairpin} or die "no hairpinRange loaded\n";
    if(($leftCenter-$hairpinRange <= $relStart)&&
       ($relStop <= $rightCenter+$hairpinRange)) {
	return 1;
    }
    return 0;
}

sub computeMaxProdHitCount {
    my($productInfo,$location,$distinctReads,$parameters) = @_;
    my $minDist = $parameters->{distanceMin} or die "failed to load parameter distanceMin\n";
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $maxCount = 0;
    my $maxProdHitCount = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    if($newType eq "miR") {
		my $sum = 0;
		my $total = 0;
		foreach my $distinctRead (@{$distinctReads->{$strand}}) {
		    # for each read...
		    my($rStart,$rStop,$readCounts,$hitCount,$libraryCounts) = @{$distinctRead};
		    # if it is associated with a product...
		    if(abs($gStart-$rStart) < $minDist) {
			# compute the average number of hits to the genome
			$sum += $hitCount * $readCounts;
			$total += $readCounts;
			#print "incrementing: ( $rStart $rStop ) $rID $hitCount->{$sample.$rID}\n";
			#print "--> $sum $total\n";
		    }
		}
		my $avg = $total ? $sum/$total : 0.0;
		#print " = $avg\n";
		if($adjProdCount > $maxCount) {
		    $maxCount = $adjProdCount;
		    $maxProdHitCount = $avg;		
		}
	    }
	 }
    }
    return $maxProdHitCount;
}

#this version of computeMaxProdHitCount does not take into account product size
sub oldComputeMaxProdHitCount {
    my($productInfo,$location,$distinctReads,$parameters) = @_;
    my $minDist = $parameters->{distanceMin} or die "failed to load parameter distanceMin\n";
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $maxCount = 0;
    my $maxProdHitCount = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	 if($newType eq "miR") {
	     my $sum = 0;
	     my $total = 0;
	     foreach my $distinctRead (@{$distinctReads->{$strand}}) {
		 # for each read...
		 my($rStart,$rStop,$readCounts,$hitCount,$libraryCounts) = @{$distinctRead};
		 # if it is associated with a product...
		 if(abs($gStart-$rStart) < $minDist) {
		     # compute the average number of hits to the genome
		     $sum += $hitCount * $readCounts;
		     $total += $readCounts;
		     #print "incrementing: ( $rStart $rStop ) $rID $hitCount->{$sample.$rID}\n";
		     #print "--> $sum $total\n";
		 }
	     }
	     my $avg = $total ? $sum/$total : 0.0;
	     #print " = $avg\n";
	     if($adjProdCount > $maxCount) {
		 $maxCount = $adjProdCount;
		 $maxProdHitCount = $avg;		
	     }
	 }
    }
    return $maxProdHitCount;
}

sub getAdjTotalReads {
    my($distinctReads,$parameters) = @_;
    my $adjTotalReads = 0;
    my $shortSize = $parameters->{"MaxShortProductSize"};
    my $longSize = $parameters->{"minLongProductSize"};
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$hitCount,$libraryCounts) = @{$dRead};
	my $length = $dStop - $dStart + 1;
	if (($length > $shortSize) && ($length < $longSize)) {
	    $adjTotalReads += $total / $hitCount;
	}
    }
    return $adjTotalReads;
}

#this version of getAdjTotalReads does not take into account product size
sub oldGetAdjTotalReads {
    my($distinctReads) = @_;
    my $adjTotalReads = 0;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$hitCount,$libraryCounts) = @{$dRead};
	$adjTotalReads += $total / $hitCount;
    }
    return $adjTotalReads;
}


sub getProductInfo {
    my($center,$fold,$basePairs,$location,$senseProducts,$parameters) = @_;
    my $test = 0;
    my ($chrom,$start,$stop,$strand) = parseLocation($location);
    my($leftCenter, $rightCenter) = @{$center};
    my @productInfo;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$center);
    foreach my $product (@{$senseProducts}) {
	my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
	my @newProductReadsData;
	my $FIRST=1;
	 my $maxOffset;
	foreach my $read (@{$productReads}) {
	    my($dRelStart,$dRelStop,$dStart,$total,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts,$rStrand) = @{$read};
	    print "dRelStart=$dRelStart\tdRelStop=$dRelStop\n" if ($test);
	    my $offset = $dRelStop - $leftCenter;
	     print "offset = $offset = $dRelStop - $leftCenter = dRelStaop - leftCenter\n" if ($test);
	    if(($offset > 0)&&($dRelStop < $rightCenter)) {
		$offset = 0;
		print "offset = $offset > 0 and drelStop = $dRelStop < $rightCenter = rightCenter\noffset set to 0" if ($test);
	    } elsif($offset > 0) {
		$offset = $dRelStart - $rightCenter;
		print "offset = $offset > 0\noffset now equals $offset = $dRelStart - $rightCenter = dRelStart - rightCenter\n" if ($test);
	    }	    
	    if (($dRelStart == $maxRelStart) && ($dRelStop == $maxRelStop)) {
		$maxOffset = $offset;
	    }
	    push(@newProductReadsData,[$dRelStart,$dRelStop,$offset,$dStart,$total,$hitCount,
				       $libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts,$rStrand]);
	}
	my $side = getProductSide($center,$leftEnd,$rightEnd,
				  $maxRelStart,$maxRelStop,$parameters);
	push(@productInfo,[$side,\@newProductReadsData,$adjProductCount,$adjMaxProductCount,
			   $maxRelStart,$maxRelStop,$maxOffset,$maxGStart,$productStrand]);	
    }
    return \@productInfo;
}

sub rebuildProducts {
    my($center,$fold,$basePairs,$productInfo,$parameters) = @_;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$center);
    my %usedMirSide;
    my %usedMorSide;
    my %newTypes;
    my @newProductInfo1;
    my @newProductInfo2;
    # sort by distance to loop, ascending order.
    my $totalReads = 0;
    foreach my $products (@{$productInfo}) {
	my($side,$productList,$adjProductCount,$adjMaxProductCount) = @{$products};
	$totalReads += $adjProductCount;
    }   
    my @sortedProductInfo = sort {abs($a->[6]) <=> abs($b->[6])} @{$productInfo};
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$adjProductCount,$adjMaxProductCount,
	   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$products};
	my $newType;
	#print "$side $relStart $relStop $offset\n";
	if($side eq "loop") {
	    $newType = "loop";
	} elsif($side eq "split") {
	    $newType = "split";
	} elsif($side eq "out") {
	    $newType = "out";
	} elsif(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if(onHairpinArm($center,$leftEnd,$rightEnd,
			    $relStart,$relStop,$parameters)) {		
		if($usedMirSide{$side}) {
		    if($usedMorSide{$side}) {
			$newType = "out";
		    } else {
			$newType = "moR";
			$usedMorSide{$side}++;
		    }
		} elsif($usedMirSide{otherSide($side)}) {
		    if(overlapsMir($center,$basePairs,$side,$relStart,$relStop,\@newProductInfo1,$parameters)) {
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
		     if($adjProductCount > 0.05*$totalReads) {
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
	} elsif(closeToHairpin($center,$relStart,$relStop,$parameters)) {
	    if($usedMirSide{$side}) {
		if($usedMorSide{$side}) {
		    $newType = "out";
		} else {
		    $newType = "moR";
		    $usedMorSide{$side}++
		}
	    } elsif($usedMirSide{otherSide($side)}) {
		if(overlapsMir($center,$basePairs,
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
			       $adjProductCount,$adjMaxProductCount,
			       $relStart,$relStop,$offset,$gStart,$productStrand]);
    }
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$adjProductCount,$adjMaxProductCount,
	   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$products};
	if(($side eq "5p")||($side eq "3p")) {
	    if($newTypes{$relStart."x".$relStop} eq "loop") {
		if(overlapsMir($center,$basePairs,$side,
			       $relStart,$relStop,\@newProductInfo1,$parameters)) {
		    $newTypes{$relStart."x".$relStop} = "miR";
		}
	    }
	}
	push(@newProductInfo2,[$side,$newTypes{$relStart."x".$relStop},$productList,
			       $adjProductCount,$adjMaxProductCount,
			       $relStart,$relStop,$offset,$gStart,$productStrand]);
    }
    return \@newProductInfo2;
}

sub addProductSize {
    my($products,$parameters) = @_;
    my $shortSize = $parameters->{"MaxShortProductSize"};
    my $longSize = $parameters->{"minLongProductSize"};
    my @reNamedProducts;
    
    foreach my $product (@{$products}) {
	my($side,$type,$productList,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $length = $relStop - $relStart + 1;
	if ($length >= $longSize) {
	    $product->[1] = "long-" . $product->[1]; #adding long to the begginging of product type
	} elsif ($length <= $shortSize) {
	    $product->[1] = "short-" . $product->[1]; #adding short to the begginging of product type
	}
	push(@reNamedProducts,$product);
    }
    return \@reNamedProducts;
}

sub overlapCurrent {
    my($productHash,$parsedRead,$parameters)=@_;
    my $minDist = $parameters->{distanceMin};
    my($readRelStart,$readRelStop,$gStart,$count,$hitCount,$libraryCounts) = @{$parsedRead};
    # sort by the first entries total count/hit count ( the most abundant entry by construction, since
    # distinctReads is already sorted by abundance. ) highest to lowest...
    my @sortedIdList = sort {$productHash->{$b}[0][3]/$productHash->{$b}[0][4] <=> $productHash->{$a}[0][3]/$productHash->{$a}[0][4]} keys(%{$productHash});
    my @idList;
    foreach my $id (@sortedIdList) {
	my($sRelStart,$sRelStop,$sGStart,$sCount,$hitCount)=@{$productHash->{$id}[0]};
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
    my($center,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;    
    my($leftCenter, $rightCenter) = @{$center};
    my $inBuffer = $parameters->{InHairpinBuffer} or die "no inHairpinBuffer loaded\n";
    my $outBuffer = $parameters->{OutHairpinBuffer} or die "no outHairpinBuffer loaded\n";   
    my $hairpinRange = $parameters->{RangeOfHairpin} or die "no hairpinRange loaded";
    if(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	if(onHairpinArm($center,$leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if($relStop <= $leftCenter+$inBuffer) {
		return "5p";
	    } elsif($rightCenter-$inBuffer<=$relStart) {	    
		return "3p";	
	    } else {
		die "unknow situation reached in getProductSide!!! ouch.\n";
	    }
	} elsif(($relStart<=$leftCenter+1)&&($rightCenter-1<=$relStop)) {
#	    print "loop named here: relStart = $relStart leftCenter = $leftCenter rightCenter = $rightCenter relStop = $relStop\n";
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

sub getAdjTotalProductReads {
    my $productInfo = shift;
    my $adjTotalReads = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) { 
	    $adjTotalReads += $adjProdCount;
	}
    }
    return $adjTotalReads;
}

#this version of getAdjTotalProductReads does not take into account product size
sub oldGetAdjTotalProductReads {
    my $productInfo = shift;
    my $adjTotalReads = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	$adjTotalReads += $adjProdCount;
    }
    return $adjTotalReads;
}

sub getMaxProductInfo {
    my $productList = shift;
    my %dReadCounts;
    my $adjProductCount = 0;
    my $adjMaxProductCount = 0;
    my($maxRelStart,$maxRelStop,$maxOffset,$maxGStart);
    foreach my $read (@{$productList}) {
        my($rRelStart,$rRelStop,$rOffset,$rGStart,$rCount,$rHitCount,$libraryCounts,$seq,$adjustedSeqCount)=@{$read};
	$dReadCounts{$rRelStart."x".$rRelStop} += $rCount * $adjustedSeqCount;
    }
    my @keys = sort {$dReadCounts{$b} <=> $dReadCounts{$a}} keys %dReadCounts;
    my $maxKey = shift(@keys);
    ($maxRelStart,$maxRelStop) = split(/x/,$maxKey); 
    foreach my $read (@{$productList}) {
        my($rRelStart,$rRelStop,$rOffset,$rGStart,$rCount,$rHitCount,$libraryCounts,$seq,$adjustedSeqCount)=@{$read};
	$adjProductCount += $rCount * $adjustedSeqCount;
	if(($rRelStart == $maxRelStart)&&($rRelStop == $maxRelStop)) {
	    $adjMaxProductCount += $rCount * $adjustedSeqCount;
	    $maxOffset = $rOffset;
	    $maxGStart = $rGStart;	    
	}
    }
    return($maxRelStart,$maxRelStop,$maxOffset,$maxGStart,$adjMaxProductCount,$adjProductCount);
}

sub getAdjTotalLibraryCounts {
    my $productList = shift;
    my %adjTotalLibraryCounts;
    foreach my $read (@{$productList}) {
	my($rRelStart,$rRelStop,$rOffset,$rGStart,$rCount,$rHitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts)=@{$read};
	foreach my $sample (keys %{$adjLibraryCounts}) {
	    $adjTotalLibraryCounts{$sample} += $adjLibraryCounts->{$sample};
	}
    }
    return \%adjTotalLibraryCounts;
}

#####################
# INPUT SUBROUTINES #
#####################

sub loadRepeatRegions {
    my $repeatRegionsFile = shift;
    my $chromSizes = shift;
    my $repeatRegions = {};
    my $fileType = checkFileFormat($repeatRegionsFile, $chromSizes);
    if ($fileType eq "fileList") {
	open(RRFL,$repeatRegionsFile) or die "could not open $repeatRegionsFile\n";
	while(<RRFL>) {
	    chomp;
	    my $repeatRegionsFile = $_;
	    my $listedFileType = checkFileFormat($repeatRegionsFile, $chromSizes);
	    $repeatRegions = readRepeatRegionsFile($repeatRegionsFile, $repeatRegions,  $listedFileType);
	}
	close(RRFL);
    } else {
	$repeatRegions = readRepeatRegionsFile($repeatRegionsFile, $repeatRegions, $fileType);
    }
    return $repeatRegions;
}

sub readRepeatRegionsFile {
    my $repeatRegionsFile = shift;
    my $repeatRegions = shift;
    my $fileType = shift;
    open(RRF, $repeatRegionsFile) or die "Could not read repeat regions file:$repeatRegionsFile\n";
    while (<RRF>) {
	if ($fileType eq 'txt' || $fileType eq 'bed') {
	    while (<RRF>) {
		chomp;
		my($chrom,$start,$stop) = split(/\t/,$_);
		push(@{$repeatRegions->{$chrom}},[$start,$stop]);
	    }
	} elsif ($fileType eq 'gff') {
	    while (<RRF>) {
		chomp;
		unless (/\#/) {
		    my ($chrom, $source, $type, $start, $stop, $score, $strand, $frame, $attribute) = split(/\t/);
		    push(@{$repeatRegions->{$chrom}},[$start,$stop]);
		}
	    }
	} else {
	    die print "unknown format of $repeatRegionsFile.  mirTRAP readRepeatRegions takes in .txt and .gff files.\n";
	}
    }
    close(RRF);
    return $repeatRegions;
}

sub readMirbaseGff3 {
    my($mirbaseGff3) = @_;
    my %hairpins;
    my %products;
    open(MBGFF3,$mirbaseGff3) or die "could not open $mirbaseGff3\n";
    while(<MBGFF3>) {
	unless(/^\#/) {
	    chomp;
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    if($type eq "miRNA_primary_transcript") {
		# hairpin region from the gff file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n"; 
		my $name = $info{Name} or die "No Name found for the line:\n$_\n"; 
		push(@{$hairpins{$chrom}},[$start,$stop,$strand,$id,$name]);
	    }
	    if($type eq "miRNA") {
		# mature product from the hairpin file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
                my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		my $parentId = $info{Derives_from} or die "No Derives_from found for the line:\n$_\n";
		push(@{$products{$parentId}},[$chrom,$start,$stop,$strand,$id,$name]);
	    }
	}
    }
    close(MBGFF3);
    return(\%hairpins,\%products);
}


sub checkFileFormat {
    my $fileName = shift;
    my $chromSizes = shift;
    my @chroms = keys %{$chromSizes};
    my $fileType;

    if (($fileType) = $fileName =~ /\.(bed|gff)$/) {
	return $fileType;
    }
    open(FP, $fileName) or die "failed to open $fileName in miRTRAP1_6::checkFileFormat";
    while (<FP>) {
	unless (/\#/) {
	    my @column = split;
	    if ($chromSizes->{$column[0]} && isInt($column[1]) && isInt($column[2])) {
		close(FP);
		return "txt";
	    } elsif (@column == 1) {
		close(FP);
		return "fileList";
	    } else {
		close(FP);
		die print "unable to determine file type for $fileName fileType\n";
	    }
	}
    }
    die print "unable to determine file type for $fileName\n";
}

sub readtRNAScanFile {
    my $tRNAScanFile = shift;
    my %tRNAs;
    my $RECORD = 0;
    if($tRNAScanFile) {
        open(TRNA,$tRNAScanFile) or die "failed to open $tRNAScanFile for reading in readtRNAScanFile()";
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
my $repeatRegionsFile = $parameters->{repeatRegionsList};
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
