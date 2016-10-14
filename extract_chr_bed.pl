#!/usr/bin/perl

=head1 NAME

extract_chr_bed.pl - Splits a standard bed file with mapped reads into smaller bed files per each chromosome 

=head1 SYNOPSIS

extract_chr_bed.pl -in all_data.bed.gz -out output_name_template -p [<pattern>|all] [-help] 

 Required arguments:
    -in       input BED or compressed BED.GZ file
    -out      output BED file name template (each newly created bed file will have following name: chrN.output_name.bed.gz)
    -p        chromosome name pattern (For example: chr1|chromosome1 e.t.c.).
    -chrNr    set the number of chromosomes to extract (default:22 - for mouse. chromosomes X, Y and M should be extracted using corresponding pattern)

 Options:
    -help -h  Help
    
 Example usage:
    extract_chr_bed.pl --input=all_data.bed --output=output_name --pattern="chrX" --chromosomes=22
	
	OR
	
    extract_chr_bed.pl -in all_data.bed -out output_name -p "chrX" -chrNr 22
	
=head1 DESCRIPTION


=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 extract_chr_bed.pl

 extract_chr_bed.pl reads standard BED file and split it by chromosomes using specified pattern. Each newly generated file for each chromosome saved using output filename template provided as one of the required parameters. Input bed file can be compressed with gzip and have *.bed.gz extension. Output file automatically saved compressed as *.bed.gz. If no gzip found files will be saved in flat text format
 
=head1 AUTHORS

=over

=item 
 Yevhen Vainshtein <yevhen.vainshtein@igb.fraunhofer.de>
 
=item 
 Vladimir Teif
 
=back

=head2 Last modified

 14 October 2016
 
=head1 LICENSE

 Copyright (C) 2012-2016 Yevhen Vainshtein, Vladimir Teif

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


use strict 'vars';
use Getopt::Long;
use Pod::Usage;
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError) ;

# Variables set in response to command line arguments
# (with defaults)

my $infile;
my $infile_name;
my $outfile="";
my $pattern;
my $needsHelp;
my $chrNr=1;

my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$infile_name,
	'output|out=s'   => \$outfile,
	'pattern|p=s' => \$pattern,
	'chromosomes|chrNr=s' => \$chrNr,
	'help|h'      => \$needsHelp
);

# Check to make sure options are specified correctly and files exist
&check_opts();

print STDERR "input BED file: $infile_name\n";
print STDERR "output BED file: $outfile\n";
print STDERR "chromosome ID template: $pattern\n";
print STDERR "total number of chromosome to extract: $chrNr\n";

my @chromosomes;
my @FHs;
my @OUT_FHs;

if ( $infile_name =~ (/.*\.gz$/) ) {
	$infile = IO::Uncompress::Gunzip->new( $infile_name )
    or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
}
else { open( $infile, "<", $infile_name ) or die "error: $infile_name cannot be opened:$!"; }

my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();
# counter for the markers we see
my $marker_count = 0;

my $regex_split_tab='.*\t(.*)';
my $regex_split_newline='\n';

my $filesize = -s $infile_name; #determine file size in bytes
my $size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "Reading $infile_name file of $filesize MBs. Please wait...\n";

my $processed_memory_size = 0;
my $offset=0;
my $not_zero_counter=0;
my $string_counter=0;
my $BUFFER_SIZE = 1024*4;

if ($chrNr > 1) {
    print STDERR "extracting $chrNr chromosomes...\n";
	
	for (my $chrID=1; $chrID <= $chrNr; $chrID++ ) {
		my $chrname = $pattern.$chrID;
		my $FH = uc($chrname);
		push (@chromosomes, $chrname);
		push (@FHs, $FH);
	}

    for (my $i=0; $i<=$#chromosomes; $i++ ) {
		push(@OUT_FHs, "*$FHs[$i]");
		my $out_filename_gz= "$chromosomes[$i]"."."."$outfile".".bed.gz";
		my $out_filename= "$chromosomes[$i]"."."."$outfile".".bed";
		$OUT_FHs[$i] = new IO::Compress::Gzip ($out_filename_gz) or open ">$out_filename" or die "Can't open $out_filename for writing: $!\n";
    }
    
    while ((my $n = read($infile, $buffer, $BUFFER_SIZE)) !=0) {
        if ($n >= $BUFFER_SIZE) {
        $buffer .= <$infile>;
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
            print STDERR "."; $processed_memory_size=0;
            }
        undef @lines;
        $buffer = "";
    }
    foreach my $OUTFH (@OUT_FHs) { close($OUTFH); } # close output files
    print STDERR "done in ", time()-$timer2, " seconds\n";
} 


else {
	my $out_filename_gz= $pattern.".".$outfile.".bed.gz";
	my $out_filename= $pattern.".".$outfile.".bed";
	my $OUT_FHs = new IO::Compress::Gzip ($out_filename_gz) or open ">$out_filename" or die "Can't open $out_filename for writing: $!\n";

    while ((my $n = read($infile, $buffer, $BUFFER_SIZE)) !=0) {
        if ($n >= $BUFFER_SIZE) {
        $buffer .= <$infile>;
        }
        my @lines = split(/$regex_split_newline/o, $buffer);
        # process each line in zone file
        foreach my $line (@lines) {
          if ($line =~ m/$pattern\s/) { print $OUT_FHs "$line\n"; }
        }
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $filesize/10) {
            print STDERR "."; $processed_memory_size=0;
	}
    undef @lines;
    $buffer = "";
    }
    
    print STDERR "done in ", time()-$timer2, " seconds\n";
    close($OUT_FHs);  
}

close($infile);


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
	if ( !-e $infile_name ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input BED file: '$infile_name!'\n"
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
