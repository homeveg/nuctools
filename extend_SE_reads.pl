#!/usr/local/bin/perl

###==================================================================================================
### Extends single-end reads by the user-defined value of the average DNA fragment length
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### extend_SE_reads.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================

use strict;
my $usage = "$0 <in.bed> <out.bed> <fragment length>\n\nInput file should be a tab-delimited file with the following columns order:\nChromosome name | start | end | col4 | col5 | strand\n";
my $in_file = shift || die $usage;
my $outfile = shift || die $usage;
my $fragment_length = shift || die $usage;

open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

my $chr_id = 0;
my $start_id = 1;
my $end_id = 2;
my $strand_id = 5;

#------------------------------------------------------------------------------
#read file with occupanicies

my $filesize = -s $in_file; #determine file size in bytes
my $size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "Reading file $in_file file of $filesize MBs. Please wait...\n";

#read file with by 4kb chanks
my $BUFFER_SIZE = 1024*4;

# open original file
open(INPUT, $in_file) or die "error: $in_file cannot be opened\n";
my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();
# counter for the markers we see
my $marker_count = 0;

my $processed_memory_size = 0;
my $offset=0;
while ((my $n = read(INPUT, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <INPUT>;
    }
    my @lines = split(/\n/, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
        my @newline=split(/\s/, $line);
        my $chr = $newline[$chr_id];
        my $start = $newline[$start_id];
        my $end = $newline[$end_id];
        my $strand = $newline[$strand_id];
        
        if ($strand =~ m/\+/) {
          $end = $start + $fragment_length;
        }
        elsif ($strand =~ m/\-/) {
          $start = $end - $fragment_length;
        }
        print OUT join("\t", $chr, $start, $end, $strand),"\n";
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

print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\ndone.\n";
close(INPUT) or die $!;
close(OUT) or die $!;
exit;
