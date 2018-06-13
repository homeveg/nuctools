#!/usr/bin/perl

=head1 NAME

 extract_rows_occup.pl - Extracts occupancy values from *.OCC.GZ file for a given genomic interval and save it as a compressed *.OCC.GZ file
 
=head1 SYNOPSIS

perl -w extract_rows_occup.pl --input=<in.bed> --output=<out.bed> --start=<selected region start> --end=<selected region end> --coordsCol=<coordinates Col.Nr.> --help [--help] 

 Required arguments:
    --input | -in      path to directory with aggregate profiles
    --output | -out    output table file name
    --start | -s       chromosomal coordinate of selected region start
    --end | -e         chromosomal coordinate of selected region end

 Options:
    --gzip | -z        compress the output
    --coordsCol | -cC  chromosomal coordinate column Nr.(default: 0)
    --help | h         Help
	
 Example usage:
 
    perl -w extract_rows_occup.pl --input=in.bed.gz --output=out.bed.gz	
	
	OR
	
    perl -w extract_rows_occup.pl -in in.bed.gz -out out.bed.gz
    
=head1 DESCRIPTION
 
=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 extract_rows_occup.pl

 extract_rows_occup.pl extracts occupancy values from *.OCC.GZ file for a given genomic interval and save it as a compressed *.OCC.GZ file
 
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

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

my $infile;
my $outfile;
my $start_interval;
my $end_interval;
my $coordsCol=0;
my $useGZ;

my $needsHelp;

my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$infile,
	'output|out=s'   => \$outfile,
	
	'start|s=s' => \$start_interval,
	'end|e=s' => \$end_interval,
	'coordsCol|cC=s'  => \$coordsCol,
	'gzip|z' => \$useGZ,

	'help|h'      => \$needsHelp
);

# Check to make sure options are specified correctly and files exist
&check_opts();

# check if GZIP is loaded
if ( ((!$ModuleGzipIsLoaded) or (!$ModuleGunzipIsLoaded)) and ($useGZ) ) {
	print STDERR "Can't work with GZIP: IO::Compress::Gzip is not on PATH\n";
	exit;
}
elsif ( (($ModuleGzipIsLoaded) and ($ModuleGunzipIsLoaded)) and ($useGZ) ) {
	print STDERR "ZGIP support enabled\n";
}
else {
	print STDERR "ZGIP support disabled\n";
	if ( ($infile =~ (/.*\.gz$/))  and (!$useGZ) ) {
		print STDERR "======================================\n";
		print STDERR "WARNING! Input file probably compressed!\n";
		print STDERR "Use --gzip parameter to enable support for file compression/decompression!";
		print STDERR "======================================\n";
		exit;
	}
}
# open pipe to Gzip or open text file for writing
  my ($gz_out_file,$out_file,$OUT_FHs);
  $out_file = $outfile;
	if ($useGZ) {
		$out_file =~ s/(.*)\.gz$/$1/;
		$gz_out_file = $out_file.".gz";
		$OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
	}
	else {
		open $OUT_FHs, '>', $outfile or die "Can't open $outfile for writing; $!\n";
	}

# open occupancy file
my $inFH;
if ( $infile =~ (/.*\.gz$/) ) {
	$inFH = IO::Uncompress::Gunzip->new( $infile )
	or die "IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip::GunzipError\n";
}
else { open( $inFH, "<", $infile ) or die "error: $infile cannot be opened:$!"; }

my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();
# counter for the markers we see
my $marker_count = 0;

my $regex_split_newline='\n';

my $filesize_in_bytes = -s $infile; #determine file size in bytes
my $size_counter_step=int($filesize_in_bytes/10);
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
        if (($newline[$coordsCol] >= $start_interval) and ($newline[$coordsCol] <= $end_interval)) {
          print $OUT_FHs join("\t",@newline), "\n";
          #print join("\t",@newline), "\n";
        }
		elsif ($newline[$coordsCol] > $end_interval) {
			last;
		}
    }
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $size_counter_step) {
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
