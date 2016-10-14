#!/usr/local/bin/perl

=head1 NAME

extract_chr_bed.pl - Splits a standard bed file with mapped reads into smaller bed files per each chromosome 

=head1 SYNOPSIS

extract_chr_bed.pl -in all_data.bed -out output_name_template -p [<pattern>|all] [-help] 

 Required arguments:
    -in       input BED file
    -out      output BED file name template (each newly created bed file will have following name: chrN.output_name.bed)
	-p        chromosome name pattern (For example: chr1|chromosome1 e.t.c.). Use "all" to extract all chromosomes from input bed

 Options:
    -help -h  Help
    
 Example usage:
    extract_chr_bed.pl -in all_data.bed -out output_name -p "chr10" [-help]
	
=head1 DESCRIPTION


=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 extract_chr_bed.pl

 extract_chr_bed.pl reads standard BED file and split it by chromosomes using specified pattern. Each newly generated file for each chromosome saved using output filename template provided as one of the required parameters 
 
=head1 AUTHORS

=over

=item 
 Yevhen Vainshtein <yevhen.vainshtein@igb.fraunhofer.de>
 
=item 
 Vladimir Teif
 
=back
 
=head1 LICENSE

 Copyright (C) 2014-2016 Yevhen Vainshtein

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

=cut

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
use Getopt::Long;
use Pod::Usage;

# Variables set in response to command line arguments
# (with defaults)

my $usage = "$0 -input=<in.bed> -output=<out.bed> -pattern=<pattern>\n";
my $infile;
my $outfile="";
my $pattern;
my $needsHelp;

my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$infile,
	'output|out=s'   => \$outfile,
	'pattern|p=s' => \$pattern,
	'help|h'      => \$needsHelp
);

# Check to make sure options are specified correctly and files exist
&check_opts();



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
	my $out_filename= "$pattern"."$outfile".".bed";
    open(OUT, ">$out_filename") || die "Can't open $out_filename for writing: $!\n";
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


# Check for problem with the options or if user requests help
sub check_opts {
	if ($needsHelp) {
		pod2usage( -verbose => 2 );
	}

	if ( !$options_okay ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Error specifying options."
		);
	}
	if ( !-e $infile ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input BED file: '$infile!'\n"
		);
	}
	if ( -e $outfile ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "'$outfile' exists in target folder\n"
		);
	}

}
