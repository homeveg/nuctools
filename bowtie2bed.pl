#!/usr/bin/perl

=head1 NAME

compare_two_conditions.pl - convert BAM/SAM or MAP formatted files to BED format

=head1 SYNOPSIS

perl -w bowtie2bed.pl.pl --input=<healthy.txt> --output=<more_than1.txt> [--verbose --help]

 Required arguments:
    --input | -i       path to input SAM/BAM/MAP file
    --output | -o      path to output BED file (OCC.GZ)
	
 Options:	
    --verbose | -v     print BED output to STDOUT
    --help | -h        Help
    
 Example usage:
    bowtie2bed.pl --input=accepte_hits.bam --output=sample.bed.gz 
	
	OR
    
    bowtie2bed.pl -i accepte_hits.bam -o sample.bed.gz 

	
=head1 DESCRIPTION

=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 bowtie2bed.pl

 bowtie2bed.pl takes as an input standard SAM, BAM or MAP file and converts to the gzip-compressed BED file. The programm require samtools installed in PATH to be able to work with BAM files
 
=head1 AUTHORS

=over

=item 
 Yevhen Vainshtein <yevhen.vainshtein@igb.fraunhofer.de>
 
=item 
 Vladimir Teif
 
=back

=head2 Last modified

 21 November 2016
 
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
use IO::Compress::Gzip qw(gzip $GzipError) ;
use File::Which;
use List::Util qw(sum);

my ($infile,$outfile,$verbose);
my $needsHelp;

my $options_okay = &Getopt::Long::GetOptions(
	'input|i=s' => \$infile,
	'output|o=s'   => \$outfile,
	
	'verbose'   => \$verbose,
	'help|h'      => \$needsHelp
);
# Check to make sure options are specified correctly and files exist
&check_opts();

# open output file (compressed id gzip installed)
my ($out_file,$gz_out_file, $OUT_FHs);
my $BAMflag;
if ( ! $verbose) {
  # if not piping to output, check if output file name specified
  if ( !$outfile )  {
  pod2usage(
    -exitval => 2,
    -verbose => 1,
    -message => "Please specify output BED file name!\n"
  );
  exit;
  }
    
  $out_file = $outfile;
  $out_file =~ s/(.*)\.gz$/$1/;
  $gz_out_file = $out_file.".gz";
  $OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
}

# open SAM/MAP file or pipe to an input BAM file
my $data;

if($infile =~ /.*\.sam$/i ) {
  print STDERR "converting SAM to BED format...\n";
  open( $data, "<", $infile ) or die "$!";
  $BAMflag="SAM";
}
elsif($infile =~ /.*\.map$/i ) {
  print STDERR "converting MAP to BED format...\n";
  open( $data, "<", $infile ) or die "$!";
  $BAMflag="MAP";
}
elsif($infile =~ /.*\.bam$/i ) {
  print STDERR "converting BAM to BED format...\n";
  my $samtools_path = which 'samtools';
  if(!$samtools_path) {
    print STDERR "can't find samtools on PATH to read from BAM file!\n";
    exit;
  }
  my $samtools_pipe = "samtools view $infile";
  open $data, "-|", $samtools_pipe or die "Pipe from $samtools_pipe failed: $!";
  $BAMflag="BAM";
}
else {  pod2usage(
    -exitval => 2,
    -verbose => 1,
    -message => "I don't know specified input format!\n"
  );
  exit;
}

while(<$data>)
{
  chomp;
	my $line = $_;
	if($line =~ /^\@.*/ ) { next; }
	
  my @read = split ' ', $line;
  my ($LENGTH, $END, $STRAND);
  my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $SEQ, $PHRED);
  
  if ( ( $BAMflag eq "BAM" ) or ( $BAMflag eq "SAM" ) ) {
    #DGKHT8Q1:204:C9G2PACXX:6:2202:13917:38047	0	chr10	3100456	0	100M	*	0	0	CACTGATAAGTGGATAATAGCCCAGAAACTTAGGATACCCAAGATATAAGATACAACTTGCCAAACGCATGAAATTCAAGAAGAACGAAGACCAAAGTGT	CCCFFFFFHGDHHIJJJJJJJIIJJJJJJJIIJIJJJJIJDIIGIJJIJJIJIJJJJJJJJJJIJJJJJIGGHHHACEFFFDEEEEDD@DDCDCBDD>C;	AS:i:1000	NM:i:0	XI:f:1	X0:i:161	XE:i:13	XR:i:100	MD:Z:100
    $QNAME = $read[0]; # read ID
    $FLAG  = $read[1]; $STRAND = $FLAG==0 ? '+' : '-' ; # strand
    $RNAME = $read[2]; # chromosome name (reference name)
    $POS   = $read[3]; # 1-based leftmost mapping POSition
    $SEQ   = $read[9]; # segment sequence
    $PHRED = $read[10];# segment Phred score
		$MAPQ = Phred($PHRED);
    $LENGTH = length($SEQ);
  }
  elsif ($BAMflag eq "MAP" ) {
    #SRR1186259.129 HWUSI-EAS270R:5:1:0:63 length=36 +       chr18   34665042        NATAATGTCACTAAGGTACACNACATACTGCAAGGA    !!!!!"""""(((#((((((!!!!!!!!!!(("(((    0       0:T>N,21:T>N
    $QNAME = $read[1]; # read ID
    $STRAND  = $read[3]; # strand
    $RNAME = $read[4]; # chromosome name (reference name)
    $POS   = $read[5]; # 1-based leftmost mapping POSition
    $SEQ   = $read[6]; # segment sequence
    $PHRED  = $read[7]; # MAPping Quality
		$MAPQ = Phred($PHRED);
    $LENGTH = length($SEQ);

  }
  
  $END = $POS + $LENGTH;

  # output BED format (https://genome.ucsc.edu/FAQ/FAQformat)
  if ($verbose) {
    print STDOUT "$RNAME\t$POS\t$END\t$QNAME\t$MAPQ\t$STRAND\n";
  }
  else {
    print $OUT_FHs "$RNAME\t$POS\t$END\t$QNAME\t$MAPQ\t$STRAND\n";
  }
}

#close($data);
if (!$verbose) {
  close($OUT_FHs);
}

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
	if ( ! -f $infile ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input file $infile: $!\n"
		);
	}

}

sub calcMean {
    return sum(@_)/@_;
}

#-------------- determine maximum value from the array ---------
sub max {
  my $max = $_[0];
  for ( @_[ 1..$#_ ] ) { $max = $_ if $_ > $max; }
  return($max);
}
#-------------- determine minimum value from the array ---------
sub min {
  my $min = $_[0];
  for ( @_[ 1..$#_ ] ) { $min = $_ if $_ < $min; }
  return($min);
}

#-------------- calculate mapping quality from Phred score ---------
sub Phred {
	  my $PHRED = $_;
    my @chars = split("", $PHRED);
    my @nums = sort { $a <=> $b } (map { unpack("C*", $_ )} @chars);
		my $MAPQ = sprintf("%.2f", calcMean(@nums)); # MAPping Quality
		return($MAPQ);
}
