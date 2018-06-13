#!/usr/bin/perl

=head1 NAME

extend_SE_reads.pl - Extends single-end reads by the user-defined value of the average DNA fragment length 

=head1 SYNOPSIS

perl -w extend_SE_reads.pl -in <in.bed> -out <out.bed> -fL <fragment length> [-cC <column Nr.> -sC <column Nr.> -eC <column Nr.> -strC <column Nr.> ] [--help] 

 Required arguments:
    --input | -in      path to directory with aggregate profiles
    --output | -out    output table file name

 Options:
 
 define column numbers in the input BED file (Nr. of the very first column is 0):
    --start_col | -sC          read start column Nr. (default: -s 1)
    --end_col | -eC            read end column Nr. (default: -e 2)
    --strand_col | -strC       strand column Nr. (default: -str 5)
    --chromosome_col | -chrC   chromosome column Nr. (default: -chr 0)
 
 additional parameters   
    --gzip | -z        compress the output
    --fragment | -fL   average DNA fragment length (default: -fL 147)
	--help | h                 Help
	
 Example usage:
 
    perl -w extend_SE_reads.pl --input=in.bed.gz --output=out.bed.gz --fragment=100 	
	
	OR
	
    perl -w extend_SE_reads.pl -in in.bed.gz -out out.bed.gz -l 100 	
    
=head1 DESCRIPTION
 
=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 average_replicates.pl

 extend_SE_reads.pl extends single-end reads by the user-defined value of the average DNA fragment length. Script works with compressed or uncompressed BED files and save output as compress *.BED.GZ

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

use strict 'vars';
use Getopt::Long;
use Pod::Usage;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

my $in_file;
my $outfile;
my $fragment_length = 147;

my $chr_id = 0;
my $start_id = 1;
my $end_id = 2;
my $strand_id = 5;

my $needsHelp;
my $useGZ;

my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$in_file,
	'output|out=s'   => \$outfile,
	
	'start_col|sC=s' => \$start_id,
	'end_col|eC=s'   => \$end_id,
	'strand_col|strC=s' => \$strand_id,
	'chromosome_col|chrC=s'   => \$chr_id,
	'fragment|fL=s'   => \$fragment_length,
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
		if ( ($in_file =~ (/.*\.gz$/))  and (!$useGZ) ) {
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

#------------------------------------------------------------------------------
#read file with occupanicies

my $filesize = -s $in_file; #determine file size in bytes
my $size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "Reading file $in_file file of $filesize MBs. Please wait...\n";

#read file with by 4kb chanks
my $BUFFER_SIZE = 1024*4;

# open occupancy file
my $inFH;
if ( $in_file =~ (/.*\.gz$/) ) {
	$inFH = IO::Uncompress::Gunzip->new( $in_file )
	or die "IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip::GunzipError\n";
}
else { open( $inFH, "<", $in_file ) or die "error: $in_file cannot be opened:$!"; }
	
my $buffer = "";
my $sz_buffer = 0;
my $timer2 = time();
# counter for the markers we see
my $marker_count = 0;

my $processed_memory_size = 0;
my $offset=0;
while ((my $n = read($inFH, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <$inFH>;
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
        print $OUT_FHs join("\t", $chr, $start, $end, $strand),"\n";
    }
    $processed_memory_size += $n;
    $offset += $n;
    if(int($processed_memory_size/1048576)>= $filesize/10) {
	print STDERR "."; $processed_memory_size=0;
	}
    undef @lines;
    $buffer = "";
}

my $duration = time()-$timer2;

print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\ndone.\n";
close($inFH) or die $!;
close($OUT_FHs) or die $!;
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
	if ( ! -e $in_file ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input BED file $in_file: $!\n"
		);
	}
	if (!$outfile ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify output BED file name\n"
		);
	}

}
