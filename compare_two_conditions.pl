#!/usr/local/bin/perl

###==================================================================================================
### Takes as input occupancy files for two conditions (which have been previously generated from 
### several replicates using average_replicates.pl), and produces two files with regions of size 
### windowSize, where the difference between the signal in condition 2 and condition 1 is 
### correspondingly larger or smaller than threshold1 and threshold 2.
###
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### compare_two_conditions.pl 
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================

use strict;
use Time::localtime;
use Time::Local;
use File::Basename;
use List::Util qw(sum);


my $usage = "$0 -input1=healthy.txt -input2=patients.txt -output1=more_than1.txt -output2=less_than1.txt -chromosome=chr1 -windowSize=100 -threshold1=0.8 -threshold2=0.5 -Col_signal=1 -Col_coord=0\n";

my ($input1,$input2);
my $Col_signal=1;
my $Col_coord=0;
my $threshold1=0.8;
my $threshold2=0.5;
my ($output1,$output2);
my $chromosome="chr1";
my $windowSize=100;
my $verbose = "no";

if (@ARGV != 0) {
    foreach my $comand_line_flag (@ARGV) {
	if ($comand_line_flag =~ /-input1=(.*)/i) { $input1 = $1; }
	if ($comand_line_flag =~ /-input2=(.*)/i) { $input2 = $1; }
	if ($comand_line_flag =~ /-output1=(.*)/i) { $output1 = $1; }
	if ($comand_line_flag =~ /-output2=(.*)/i) { $output2 = $1; }
	if ($comand_line_flag =~ /-Col_signal=(.*)/i) { $Col_signal = $1; }
	if ($comand_line_flag =~ /-Col_coord=(.*)/i) { $Col_coord = $1; }
	if ($comand_line_flag =~ /-threshold1=(.*)/i) { $threshold1 = $1; }
	if ($comand_line_flag =~ /-threshold2=(.*)/i) { $threshold2 = $1; }
	if ($comand_line_flag =~ /-chromosome=(.*)/i) { $chromosome = $1; }
	if ($comand_line_flag =~ /-windowSize=(.*)/i) { $windowSize = $1; }
	if ($comand_line_flag =~ /--verbose/i) { $verbose = "yes"; }
    }
}
else {
    print STDERR $usage,"\n";
    exit;
}

print STDERR "======================================\n";
print STDERR "input file 1:",$input1, "\n";
print STDERR "input file 2:",$input2, "\n";
print STDERR "output file 1:",$output1, "\n";
print STDERR "output file 2:",$output2, "\n";
print STDERR "======================================\n";
print STDERR "Occupancy column ID: ",$Col_signal, "\n";
print STDERR "Coordinates column ID: ",$Col_coord, "\n";
print STDERR "Upper threshold: ",$threshold1, "\n";
print STDERR "Lower threshold: ",$threshold2, "\n";
print STDERR "======================================\n";
print STDERR "chromosome name: ",$chromosome, "\n";
print STDERR "bin size: ",$windowSize, "\n";
print STDERR "======================================\n";
print STDERR "print all data to STDOUT: ",$verbose, "\n";

my $tm = localtime;
print STDERR "-----------------------\n",join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),"\n-----------------------\n";

my $filesize = -s $input1; #determine file size in bytes
my $size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "Reading $input1 file of $filesize MBs. Please wait...\n";

#read file with by 4kb chanks
#@coord_occ_array=();
my $BUFFER_SIZE = 1024*4;

# open original file
open(INPUT, $input1) or die "error: $input1 cannot be opened\n";
my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();

my $regex_split_newline='\n';

my $processed_memory_size = 0;
my $offset=0;

my (@occup1,@pos1);

while ((my $n = read(INPUT, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <INPUT>;
    }
    my @lines = split(/$regex_split_newline/o, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
        my @string;
        push (@string, split("\t",$line));
        push (@occup1,$string[$Col_signal]);
        push (@pos1,$string[$Col_coord]);
        undef @string;
    }
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $filesize/10) {
        print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
        }
    undef @lines;
    $buffer = "";
}

my $duration = time()-$timer2;

print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\ndone.\n\n";
close(INPUT) or die $!;



$filesize = -s $input2; #determine file size in bytes
$size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "Reading $input2 file of $filesize MBs. Please wait...\n";

# open original file
open(INPUT, $input2) or die "error: $input2 cannot be opened\n";
$buffer = "";
$sz_buffer = 0;
$timer2 = time();

$processed_memory_size = 0;
$offset=0;

my (@occup2,@pos2);

while ((my $n = read(INPUT, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <INPUT>;
    }
    my @lines = split(/$regex_split_newline/o, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
        my @string;
        push (@string, split("\t",$line));
        push (@occup2,$string[$Col_signal]);
        push (@pos2,$string[$Col_coord]);
        undef @string;
    }
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $filesize/10) {
        print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
        }
    undef @lines;
    $buffer = "";
}

$duration = time()-$timer2;

print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\ndone.\n";
close(INPUT) or die $!;

open(OUT1,">$output1") or die $!;
open(OUT2,">$output2") or die $!;
print STDERR "\n======================\nstart filtering...";
my $above_counter=0;
my $below_counter=0;
my $null_counter=0;

for (my $i=1; $i<=$#pos1; $i++) {
    my $norm_difference;
    if (($occup1[$i]==0) or ($occup2[$i])==0) { $null_counter++; next; }

    if (($occup1[$i]+$occup2[$i])==0) { $norm_difference=0; }
    else { $norm_difference=2*($occup1[$i]-$occup2[$i])/($occup1[$i]+$occup2[$i]); }
    my $start_region = $pos1[$i]-$windowSize;
    my $end_region = $pos1[$i];
    my $above_below_flag="between";
    if ($norm_difference > $threshold1) {
        print OUT1 join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1[$i], $occup2[$i]),"\n";
	$above_counter++;
	$above_below_flag="above";
    }
    if ($norm_difference < $threshold2) {
        print OUT2 join("\t",$chromosome, $start_region, $end_region, $norm_difference, $occup1[$i], $occup2[$i]),"\n";
	$below_counter++;
	$above_below_flag="below";
    }
    
    if ($verbose eq "yes") {
	#code
	print STDERR join("\t",$above_below_flag, $chromosome, $start_region, $end_region, $norm_difference, $occup1[$i],$occup2[$i]), "\n";
    }
    
}

my $between=$#pos1-$above_counter-$below_counter-$null_counter;
print STDERR "done!\n======================\n";
print STDERR "Upper threshold: ",$threshold1, "\n";
print STDERR "Lower threshold: ",$threshold2, "\n";
print STDERR "======================\n",
	     "$null_counter from $#pos1 entries eq to 0 in both datasets\n",
	     "$between from $#pos1 entries are between 2 thresholds\n",
	     "$above_counter from $#pos1 entries saved to $output1\n",
	     "$below_counter from $#pos1 entries saved to $output2\n======================\n";
close(OUT1) or die $!;
close(OUT2) or die $!;
exit;
