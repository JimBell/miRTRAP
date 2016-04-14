#!/usr/bin/perl -w
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 <sam file>\n";

my $samFile = $ARGV[0] or die $USAGE;
my($fileBase) = $samFile =~ /([^\/]*)\.sam/;
my $outputFastaFile = $fileBase . "_condensed.fasta";

createFasta($samFile,$outputFastaFile);

sub createFasta {
    my($samFile,$outputFasta) = @_;
#    my @samArgs = ("cat",$samFile,"|","awk","'!/^@/'","|","awk","'{print \"\$1\t\$10\"}'","|","sort","|","uniq","|","awk","'{print \$2}'","|","sort","|","uniq","-c");
#    open(SAM,"-|",@samArgs);
    my $command = "cat $samFile | awk '!\(/^@/\) {print \"\$1\t\$10\"}' > $samFile.temp";
    system("cat $samFile | awk '!\(/^@/\) {print \"".'$1'."\t".'$10'."\"}' > $samFile.temp");# | sort | uniq | awk '{print \$2}' | sort | uniq -c > $samFile.temp");
#    open(OFSTA,">$outputFasta");
#    my $readNum = 1;
#    while (<SAM>) {
#	chomp;
#	print $_."\n";
#	my($count,$seq) = split(/\t/);
#	print OFSTA ">dr".$readNum."_x".$count."\n";;
#	print OFSTA "$seq\n";
#	$readNum++;
 #   }
#    close(SAM);
#    close(OFSTA);
}
