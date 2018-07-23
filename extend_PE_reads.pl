#!/usr/bin/perl

=head1 NAME

extend_PE_reads.pl - Takes as input bed file with mapped paired-end reads (two lines per paired read) and reformat it by creating a smaller bed file
with one line per nucleosome in the following format:
(1) chromosome, (2) nucleosome start, (3) nucleosome end, (4) nucleosome length

=head1 SYNOPSIS

perl -w extend_PE_reads.pl -in <in.bed> -out <out.bed> [--verbose --help] 

 Required arguments:
    --input | -in      path to input bed file with mapped paired-end reads
    --output | -out    output "one-line-per-paired-end-reads" bed file name

 Options:
    --NucLength | -nl  maximum expected read length (default: 1000)
    --gzip | -z        compress the output
    --verbose | -v     print converted output to STDOUT
    --fragment | -fL   average DNA fragment length (default: -fL 147)
	--extend | -e      extend fragments to the length defined by a "--fragment" option (default: -fL 147)

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
use List::Util qw(min max);
use Time::localtime;
use Time::Local;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

my $infile;
my $outfile;

my $needsHelp;
my $useGZ;
my $verbose;
my $NucLength = 1000;
my $fragment_length = 147;
my $extendFragment;

my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$infile,
	'output|out=s'   => \$outfile,
    'NucLength|nl=s' => \$NucLength,
	'gzip|z' => \$useGZ,
	'verbose'   => \$verbose,
	'fragment|fL=s'   => \$fragment_length,
	'extend|e'      => \$extendFragment,

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



#  Time count Initialization
my $timer1=time();
my $tm = localtime;
my $start_sec = $tm -> [0];
my $start_min = $tm ->[1];
my $start_hour = $tm ->[2];
my $start_time = time();

# Display input parameters
print STDERR "======================================\n";
print STDERR "Started:\t$start_hour:$start_min:$start_sec\n";
print STDERR "======================================\n";
print STDERR "in file:",$infile, "\n";
print STDERR "out file:",$outfile, "\n";
print STDERR "maximum fragment length: ",$NucLength, "\n";
if ( defined $extendFragment) {
	print STDERR "extend all short fragments symmetrically to an expected fragment length: $fragment_length\n";
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
while ((my $n = read($inFH, $buffer, $BUFFER_SIZE)) !=0) {
    if (($n >= $BUFFER_SIZE) or (($n == $filesize_in_bytes))) {
        $buffer .= <$inFH>;
    }
    my @lines = split(/\n/, $buffer);    
    my $end_index=$#lines;
    for (my $i=0; $i<=$end_index; $i++) {
      	my ($line1,$line2);
	if($last_line) {
	  unshift(@lines, $last_line);
	  $end_index=$#lines;
	  undef $last_line;
	}
	if(($i==$end_index) && ($end_index % 2 == 0) && ($lines[$#lines] =~ /^chr.*/ ))  { $last_line= $lines[$#lines]; last; }
	$line1=$lines[$i]; chomp($line1);
	$line2=$lines[$i+1]; chomp($line2);
	
	my @newline1=split(/\t/, $line1);
	my @newline2=split(/\t/, $line2);
	
	my $chr_name_1=$newline1[0];
	my $chr_name_2=$newline2[0];
    
    my $min = min ($newline1[1],$newline2[1]);
    my $max = max ($newline1[2],$newline2[2]);
    my $nuc_length = $max - $min;

    my $read_1=$newline1[3];
	my $read_2=$newline2[3];
    
   	if (($read_1 eq $read_2) & ($chr_name_1 eq $chr_name_2) & ($nuc_length >0) & ($nuc_length < $NucLength))  {
		if ( defined $extendFragment & ($nuc_length<$fragment_length) ){
			my $delta=int ($fragment_length-$nuc_length)/2;
			$min-=$delta;$max+=$delta;$nuc_length=$max - $min;
		}
        print $OUT_FHs join("\t", $chr_name_1, $min, $max, $nuc_length), "\n";
        if ($verbose) {
            print STDOUT join("\t", $chr_name_1, $min, $max, $nuc_length), "\n";
        }
	}

	$i++;
    }
    if($#lines % 2)  {
      undef $last_line;
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
print STDERR "job finished! Bye!\n";
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

}
