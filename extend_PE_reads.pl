#!/usr/local/bin/perl

###==================================================================================================
### Takes as input bed file with mapped paired-end reads (two lines per paired read) and reformats 
### it by creating a smaller bed file with one line per nucleosome in the following format: 
### (1) chromosome, (2) nucleosome start, (3) nucleosome end, (4) nucleosome length. 
###
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### extend_PE_reads.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================

use strict;
my $usage = "$0 <in.bed> <out.bed>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
print "script started... \n";

open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";
my @stuff=();


my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();
# counter for the markers we see
my $marker_count = 0;

my $filesize_in_bytes = -s $infile; #determine file size in bytes
my $size_counter_step=int($filesize_in_bytes/100);
my $filesize = int($filesize_in_bytes/1048576); # filesize in megabytes

print STDERR "Reading $infile file of $filesize MBs. Please wait...\n";
my $processed_memory_size = 0;
my $offset=0;
my $not_zero_counter=0;
my $string_counter=0;
my $BUFFER_SIZE = 1024;
my $old_coordinate=1;
my $last_line;
while ((my $n = read(IN, $buffer, $BUFFER_SIZE)) !=0) {
    if (($n >= $BUFFER_SIZE) or (($n == $filesize_in_bytes))) {
        $buffer .= <IN>;
    }
    my @lines = split(/\n/, $buffer);
    # process each line in zone file
    
    my $end_index=$#lines;
    for (my $i=0; $i<=$end_index; $i++) {
      	my ($line1,$line2);
	if($last_line) {
	  unshift(@lines, $last_line);
	  $end_index=$#lines;
	  undef $last_line;
	}
	if(($i==$end_index) && ($end_index % 2 == 0) && ($lines[$#lines] =~ /^chr.*/ )) 
	    { $last_line= $lines[$#lines]; last; }
	$line1=$lines[$i]; chomp($line1);
	$line2=$lines[$i+1]; chomp($line2);
	
	my @newline1=split(/\t/, $line1);
	      #print "newline1:",$#newline1," ",join("\t", @newline1), "\n";
	my @newline2=split(/\t/, $line2);
	     #print "newline2:",$#newline2," ", join("\t", @newline2), "\n";
	
	
	my $chr_name=$newline1[0];
	my $start_nuc=$newline1[1];
	my $end_nuc=$newline2[2];
	my $nuc_length=$end_nuc-$start_nuc;
	
	 #
	#print join("\t", $chr_name, $start_nuc, $end_nuc, $nuc_length), "\n";
	print OUT join("\t", $chr_name, $start_nuc, $end_nuc, $nuc_length), "\n";
	
	$i++;
    }
    if($#lines % 2)  {
      undef $last_line;
    }
    
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $filesize/100) {
        print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
        #last;
        }
    undef @lines;
    $buffer = "";
}
close(IN);
close(OUT);
print STDERR "job finished! Bye!\n";
exit;