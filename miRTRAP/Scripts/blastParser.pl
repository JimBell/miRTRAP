#!/usr/local/bin/perl -w
use Bio::SearchIO;
$|=1;
my($blastFile,@params)=@ARGV;
my $eValueThreshold=100;
my $scoreThreshold=0;
my $fracConservedThresh=0;
my $SAMEORI=0;
my $PRINTTOP=0;

readParams(@params);
parseBlastFile($blastFile);

sub readParams {
    my(@params)=@_;
    my($i,$flag,$rest,@temp,$number,@a);
    @a=split(/-/,join(" ",@params));
    for($i=1;$i<@a;$i++) {
        if($a[$i]) {
            $flag=substr($a[$i],0,1);
            $rest=substr($a[$i],1,length($a[$i])-1);
	    if($flag eq "s") {
		$rest =~ tr/\_/\-/;
                $scoreThreshold=$rest;
            } elsif($flag eq "e") {
		$rest =~ tr/\_/\-/;
                $eValueThreshold=$rest;
	    } elsif($flag eq "f") {
                $fracConservedThresh=$rest;
            } elsif($flag eq "o") {
		$SAMEORI=1;
	    } elsif($flag eq "t") {
		$PRINTTOP=1;
	    }
        }
    }
}

sub parseBlastFile {
    my($blastFile)=@_;
    my($maxScore,$maxHit,$maxHSP);
    my $input = new Bio::SearchIO(-format => 'blast',
				  -file => $blastFile);
    while(my $query = $input->next_result) {
	$maxScore=0;
	while(my $hit = $query->next_hit) {
	    while(my $hsp = $hit->next_hsp) {
		if(($hsp->evalue <= $eValueThreshold)&&($hsp->score >= $scoreThreshold)&&($hsp->frac_conserved >= $fracConservedThresh)) {
		    unless(($SAMEORI)&&($hsp->query->strand eq $hsp->hit->strand)) {
			if($PRINTTOP) {
			    if($hsp->score > $maxScore) {
				$maxScore=$hsp->score;
				$maxHit=$hit;
				$maxHSP=$hsp;
			    }
			} else {
			    print $query->query_name,":",  $hsp->query->start,"-", $hsp->query->end,":", blastOriMap($hsp->query->strand)," hits "; 
			    print $hit->name,":", $hsp->hit->start,"-", $hsp->hit->end,":", blastOriMap($hsp->hit->strand)," "; 
			    print $hsp->evalue," ", $hsp->score," ", $hsp->length, " ",$hsp->frac_conserved,"\n";   
			}
		    }
		}
	    }
	}
	if($PRINTTOP) {
	    print $query->query_name,":",  $maxHSP->query->start,"-", $maxHSP->query->end,":", blastOriMap($maxHSP->query->strand)," hits "; 
	    print $maxHit->name,":", $maxHSP->hit->start,"-", $maxHSP->hit->end,":", blastOriMap($maxHSP->hit->strand)," "; 
	    print $maxHSP->evalue," ", $maxHSP->score," ", $maxHSP->length, " ",$maxHSP->frac_conserved,"\n";    
	}
    }
}

sub blastOriMap {
  my($ori)=@_;
  if($ori eq "-1") {
    return "-";
  } elsif($ori eq "1") {
    return "+";
  }
}
