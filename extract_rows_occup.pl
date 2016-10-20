#!/usr/bin/perl

=head1 NAME

 extract_rows_occup.pl - Extracts occupancy values from *.OCC.GZ file for a given genomic interval and save it as a compressed *.OCC.GZ file
 
=head1 SYNOPSIS

perl -w extract_rows_occup.pl -input=<in.bed> -output=<out.bed> -start=<selected region start> -stop=<selected region end> --help [--help] 

 Required arguments:
    --input | -in      path to directory with aggregate profiles
    --output | -out    output table file name
    --start | -s       chromosomal coordinate of selected region start
    --end | -e         chromosomal coordinate of selected region end

 Options:
	--help | h                 Help
	
 Example usage:
 
    perl -w extend_PE_reads.pl --input=in.bed.gz --output=out.bed.gz	
	
	OR
	
    perl -w extend_PE_reads.pl -in in.bed.gz -out out.bed.gz
    
=head1 DESCRIPTION
 
=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 average_replicates.pl

 extend_PE_reads.pl Takes as input bed file with mapped paired-end reads (two lines per paired read) and reformat it by creating a smaller bed file with one line per nucleosome in the following format:
(1) chromosome, (2) nucleosome start, (3) nucleosome end, (4) nucleosome length

=head1 AUTHORS

=over

=item 
 Yevhen Vainshtein <yevhen.vainshtein@igb.fraunhofer.de>
 
=item 
 Vladimir Teif
 
=back

=head2 Last modified

 18 October 2016
 
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

use strict;
use Getopt::Long;
use Pod::Usage;
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);

my $infile;
my $outfile;
my $start_interval;
my $end_interval;

my $needsHelp;

my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$infile,
	'output|out=s'   => \$outfile,
	
	'start|p=s' => \$start_interval,
	'end|cC=e' => \$end_interval,

	'help|h'      => \$needsHelp
);

# Check to make sure options are specified correctly and files exist
&check_opts();

# open pipe to Gzip or open text file for writing
my $out_file = $outfile;
$out_file =~ s/(.*)\.gz$/$1/;
my $gz_out_file = $out_file.".gz";
my $OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";

# open occupancy file
my $inFH;
if ( $infile =~ (/.*\.gz$/) ) {
	$inFH = IO::Uncompress::Gunzip->new( $infile )
	or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
}
else { open( $inFH, "<", $infile ) or die "error: $infile cannot be opened:$!"; }

my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();
# counter for the markers we see
my $marker_count = 0;

my $regex_split_newline='\n';

my $filesize_in_bytes = -s $infile; #determine file size in bytes
my $size_counter_step=int($filesize_in_bytes/100);
my $filesize = int($filesize_in_bytes/1048576); # filesize in megabytes

print STDERR "Reading $infile file of $filesize MBs. Please wait...\n";
my $processed_memory_size = 0;
my $offset=0;
my $BUFFER_SIZE = 1024;
my $old_coordinate=1;

while ((my $n = read($inFH, $buffer, $BUFFER_SIZE)) !=0) {
    if (($n >= $BUFFER_SIZE) or (($n == $filesize_in_bytes))) {
        $buffer .= <$inFH>;
    }
    my @lines = split(/$regex_split_newline/o, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
        chomp($line);
        my @newline=split(/\s/, $line);
        if (($newline[0] >= $start_interval) and ($newline[0] <= $end_interval)) {
          print $OUT_FHs join("\t",@newline), "\n";
          print join("\t",@newline), "\n";
        }
		elsif ($newline[0] > $end_interval) {
			last;
		}
    }
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $filesize/100) {
        print STDERR "."; $processed_memory_size=0;
        }
    undef @lines;
    $buffer = "";
}
close($inFH);
close($OUT_FHs);

print STDERR "done!\nJob finished in ", time()-$timer2, " seconds.\n";
exit;

#--------------------------------------------------------------------------
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
	if ( ! -e $infile ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input BED file $infile: $!\n"
		);
	}
	if (!$outfile ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify output BED file name\n"
		);
	}
	if (!$start_interval ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify starting coordinate of the interval\n"
		);
	}
	if (!$end_interval ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify end coordinate of the interval\n"
		);
	}

}
