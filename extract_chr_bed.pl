#!/usr/local/bin/perl

###==================================================================================================
### Splits a standard bed file with mapped reads into smaller bed files per each chromosome
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### extract_chr_bed.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================

use strict 'vars';
my $usage = "$0 -input=<in.bed> -output=<out.bed> -pattern=<pattern>\n";
my $infile;
my $outfile="";
my $pattern;

if (@ARGV != 0) {
    foreach my $comand_line_flag (@ARGV) {
	if ($comand_line_flag =~ /-input=(.*)/i) { $infile = $1; }
	if ($comand_line_flag =~ /-output=(.*)/i) { $outfile = $1;}
	if ($comand_line_flag =~ /-pattern=(.*)/i) { $pattern = $1;}
    }
}
else { die $usage; }



print STDERR "infile: $infile\n";
print STDERR "outfile: $outfile\n";
print STDERR "chromosome: $pattern\n";

my @chromosomes = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY", "chrM");
my @FHs = ("CHR1","CHR2","CHR3","CHR4","CHR5","CHR6","CHR7","CHR8","CHR9","CHR10","CHR11","CHR12","CHR13","CHR14","CHR15","CHR16","CHR17","CHR18","CHR19","CHR20","CHR21","CHR22","CHRX","CHRY", "CHRM");
my @OUT_FHs;

open(IN, $infile) or die "error: $infile cannot be opened: $!\n";
my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();
# counter for the markers we see
my $marker_count = 0;

my $regex_split_tab='.*\t(.*)';
my $regex_split_newline='\n';

my $filesize = -s $infile; #determine file size in bytes
my $size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "Reading $infile file of $filesize MBs. Please wait...\n";

my $processed_memory_size = 0;
my $offset=0;
my $not_zero_counter=0;
my $string_counter=0;
my $BUFFER_SIZE = 1024*4;
my @origin; # to keep source data



if ($pattern eq "all") {
    print STDERR "extracting all chromosomes...\n";
    for (my $i=0; $i<=$#chromosomes; $i++ ) {
      my $out_filename= "$chromosomes[$i]"."$outfile".".bed";
      push(@OUT_FHs, "*$FHs[$i]");
      open($OUT_FHs[$i], ">$out_filename") || die "Can't open $out_filename for writing: $!\n";
    }
    
    while ((my $n = read(IN, $buffer, $BUFFER_SIZE)) !=0) {
        if ($n >= $BUFFER_SIZE) {
        $buffer .= <IN>;
        }
        my @lines = split(/$regex_split_newline/o, $buffer);
        # process each line in zone file
        foreach my $line (@lines) {
            for (my $i=0; $i<=$#chromosomes; $i++ ) {
                my $chromosome_id = $chromosomes[$i];
                my $OUTFH = $OUT_FHs[$i];
                if ($line =~ m/$chromosome_id\s/) { print $OUTFH "$line\n"; last; }
            }
        }
        $processed_memory_size += $n;
        $offset += $n;
        if(int($processed_memory_size/1048576)>= $filesize/10) {
            print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
            #last;
            }
        undef @lines;
        $buffer = "";
    }
    foreach my $OUTFH (@OUT_FHs) { close($OUTFH); } # close output files
    print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\nDone!\n";
} 


else {
    open(OUT, ">$outfile") || die "Can't open $outfile for writing: $!\n";
    while ((my $n = read(IN, $buffer, $BUFFER_SIZE)) !=0) {
        if ($n >= $BUFFER_SIZE) {
        $buffer .= <IN>;
        }
        my @lines = split(/$regex_split_newline/o, $buffer);
        # process each line in zone file
        foreach my $line (@lines) {
          if ($line =~ m/$pattern\s/) { print OUT "$line\n"; }
        }
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $filesize/10) {
        print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
        #last;
        }
    undef @lines;
    $buffer = "";
    }
    
    print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\nDone!\n";
    close(OUT);  
}

close(IN);
