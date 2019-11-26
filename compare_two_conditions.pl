#!/usr/bin/perl

=head1 NAME

compare_two_conditions.pl - identify regions with highest variance between control and sample condition (based on replicates) 

=head1 SYNOPSIS

perl -w compare_two_conditions.pl --input1=<control.occ> --input2=<experimental.occ> --output1=<more_than1.txt> --output2=<less_than1.txt> --chromosome="chr1" [--windowSize=100 --threshold1=0.8 --threshold2=0.5 --allowNull --Col_signal=1 --Col_coord=0 --Col_StDev=2 --Col_RelErr=3 --withError --gzip --verbose --help]

 Required arguments:
    --input1 | -i1        input average occupancy file or OCC.GZ extended file, condition I (the output file of average_replicates.pl script)
    --input2 | -i2        input average occupancy file or OCC.GZ extended file, condition II (the output file of average_replicates.pl script)
    --output1 | -o1       output high/low varience regions file (OCC.GZ)
    --output2 | -o2       output high/low varience regions file (OCC.GZ)
	
 Options:
    define column numbers in the input occupancy files (Nr. of the very first column is 0):
    --Col_coord | -cC           coordinate column Nr. (default: 0)
    --Col_signal | -sC          occupancy column Nr. (default: 1)
    --Col_StDev | -cS           Standard Deviation column Nr. (default: 2)
    --Col_RelErr | -cR          Relative Error column Nr. (default: 3)
	
   additional parameters
    --chromosome | -c           chromosome ID (mandatory parameter)
    --windowSize | -w           running window size. Use same value as for average occupancy calculation (default: 1)
    --threshold1 | -t1          upper threshold (default: 0.95)
    --threshold2 | -t2          lower threshold (default: -0.95)
	--allowNull | -aN           allow occupancy 0 in either conditions
	--withError | -wE           print StDev and RelError to the output

	
    --verbose | -v              verbose output 
    --gzip | -z                 compress the output
    --help | -h                 Help
    
 Example usage:
    compare_two_conditions.pl --input1=healthy.occ --input2=patients.occ --output1=more_than0.8.txt --output2=less_than0.5.txt --window=1000 --threshold1=0.5 --threshold2=0.8
	
	OR
    
    compare_two_conditions.pl -i1 healthy.occ -i2 patients.occ -o1 more_than0.8.txt -o2 less_than0.5.txt -w 1000 -t1 0.5 -t2 0.8

 Note:
    default column numbers refer to a standard bed file format:
	
	chromosome | region start | region end | read ID | score | strand
	=============================================================
	
	output file contains following columns (default):
	
	chromosome | region start | region end | normalized differnece | occupancy cond. 1 | occupancy cond. 2
	======================================================================================================
	
	output file contains following columns (with the --withError flag active):
	
	chromosome | region start | region end | normalized differnece | occupancy cond. 1 | standard deviation cond. 1 | rel.error cond. 1 | occupancy cond. 2 | standard deviation cond. 2 | rel.error cond. 2 
	======================================================================================================

	
=head1 DESCRIPTION

=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 compare_two_conditions.pl

 compare_two_conditions.pl takes as an input average occupancy files for two conditions (which have been previously generated from several replicates using average_replicates.pl), and produces two files with regions of size windowSize, where the difference between the signal in condition 2 and condition 1 is correspondingly larger or smaller than threshold 1 and threshold 2. The normalized difference value between average ocuupanies could be in the range from -1 to +1. Values below lower threshold indicating decrease of nucleosom occupancy in corresponding region in sample 1 comparing to sample 2. Values above the threshold indicating the increase of nucleosome occupancy.  

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
use Time::localtime;
use Time::Local;
use File::Basename;
use List::Util qw(sum);
use POSIX;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

my ($input1,$input2);
my $Col_signal=1;
my $Col_coord=0;
my $Col_StDev=2;
my $Col_RelErr=3;
my $threshold1=0.95;
my $threshold2=-0.95;
my ($output1,$output2);
my $chromosome;
my $windowSize=1;
my $verbose;
my $needsHelp;
my $useGZ;
my $allowNull;
my $withError;


my $options_okay = &Getopt::Long::GetOptions(
	'input1|i1=s' => \$input1,
	'input2|i2=s' => \$input2,
	'output1|o1=s'   => \$output1,
	'output2|o2=s'   => \$output2,
	
	'chromosome|c=s'  => \$chromosome,
	'windowSize|w=s' => \$windowSize,
	
	'threshold1|t1=s' => \$threshold1,
	'threshold2|t2=s' => \$threshold2,
	'allowNull|aN' => \$allowNull,
	
	'Col_signal|sC=s' => \$Col_signal,
	'Col_coord|cC=s' => \$Col_coord,
	'Col_StDev|cS=s' => \$Col_StDev,
	'Col_RelErr|cR=s' => \$Col_RelErr,
	'withError|wE'    => \$withError,
	'verbose|v'   => \$verbose,
	'gzip|z' => \$useGZ,
	
	'help|h'      => \$needsHelp
);

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
	if ( ($input1 =~ (/.*\.gz$/))  and (!$useGZ) ) {
		print STDERR "======================================\n";
		print STDERR "WARNING! Input file probably compressed!\n";
		print STDERR "Use --gzip parameter to enable support for file compression/decompression!";
		print STDERR "======================================\n";
		exit;
	}
}

print STDERR "======================================\n";
print STDERR "input file 1:",$input1, "\n";
print STDERR "input file 2:",$input2, "\n";
print STDERR "output file 1:",$output1, "\n";
print STDERR "output file 2:",$output2, "\n";
print STDERR "======================================\n";
print STDERR "Occupancy column ID: ",$Col_signal, "\n";
print STDERR "Coordinates column ID: ",$Col_coord, "\n";
print STDERR "Upper threshold: ",$threshold1, "\n";
print STDERR "Lower threshold: ",$threshold2, "\n";
if($allowNull) {
print STDERR "Allow zero in one of the conditions\n";
}
print STDERR "======================================\n";
if($chromosome) {
print STDERR "chromosome name: ",$chromosome, "\n";
} else {
	print STDERR "please specify chromosome ID. Exiting...\n";
	exit 1;
}
print STDERR "bin size: ",$windowSize, "\n";
print STDERR "======================================\n";
print STDERR "print all data to STDOUT: ",$verbose, "\n";

my $tm = localtime;
print STDERR "-----------------------\n",join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),"\n-----------------------\n";


my %occupancy;

ReadFile($input1, 1, 2, $Col_coord, $Col_signal,$Col_StDev,$Col_RelErr, \%occupancy);
ReadFile($input2, 2, 1, $Col_coord, $Col_signal,$Col_StDev,$Col_RelErr, \%occupancy);

# open pipe to Gzip or open text file for writing
my ($out_file,$gz_out_file, $OUT1_FHs, $OUT2_FHs);
$out_file = $output1;
if ($useGZ) {
	$out_file =~ s/(.*)\.gz$/$1/;
	$gz_out_file = $out_file.".gz";
	$OUT1_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
}
else {
	open $OUT1_FHs, '>', $output1 or die "Can't open $output1 for writing; $!\n";
}

$out_file = $output2;
if ($useGZ) {
	$out_file =~ s/(.*)\.gz$/$1/;
	$gz_out_file = $out_file.".gz";
	$OUT2_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
}
else {
	open $OUT2_FHs, '>', $output2 or die "Can't open $output2 for writing; $!\n";
}

print STDERR "\n======================\nstart filtering...";
my $above_counter=0;
my $below_counter=0;
my $null_counter=0;

for my $position ( sort {$a<=>$b} keys %occupancy) {
	my $occup1 = $occupancy{$position}{1};
	my $occup2 = $occupancy{$position}{2};
	
	my $stdev1 = $occupancy{$position}{1}{'stdev'};
	my $stdev2 = $occupancy{$position}{2}{'stdev'};

	my $relerr1 = $occupancy{$position}{1}{'relerr'};
	my $relerr2 = $occupancy{$position}{2}{'relerr'};
	
	my $norm_difference= ($occup1-$occup2)/($occup1+$occup2);

    my $start_region = $position-$windowSize;
    my $end_region = $position;
    my $above_below_flag="between";
    if (($norm_difference > $threshold1) &&  ($occup1 > 0) && ($occup2 > 0)) {
		if($withError) {
			print $OUT1_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $stdev1, $relerr1, $occup2, $stdev2, $relerr2),"\n";
		} else {
			print $OUT1_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $occup2),"\n";
		}
		$above_counter++;
		$above_below_flag="above";
    }
    elsif (($norm_difference > $threshold1) &&  ($allowNull)) {
		if($withError) {
			print $OUT1_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $stdev1, $relerr1, $occup2, $stdev2, $relerr2),"\n";
		} else {
			print $OUT1_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $occup2),"\n";
		}
		$above_counter++;
		$above_below_flag="above";
    }
	
    if (($norm_difference < $threshold2) &&  ($occup1 > 0) && ($occup2 > 0)){
		if($withError) {
			print $OUT2_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $stdev1, $relerr1, $occup2, $stdev2, $relerr2),"\n";
		} else {
			print $OUT2_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $occup2),"\n";
		}
		$below_counter++;
		$above_below_flag="below";
    }
    elsif (($norm_difference < $threshold2) &&  ($allowNull)){
		if($withError) {
			print $OUT2_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $stdev1, $relerr1, $occup2, $stdev2, $relerr2),"\n";
		} else {
			print $OUT2_FHs join("\t",$chromosome, $start_region , $end_region, $norm_difference , $occup1, $occup2),"\n";
		}
		$below_counter++;
		$above_below_flag="below";
    }
    if ($verbose) {
	#code
		if($withError) {
			print STDERR join("\t",$above_below_flag, $chromosome, $start_region, $end_region, $norm_difference, $occup1,$stdev1, $relerr1, $occup2, $stdev2, $relerr2), "\n";
		} else {
			print STDERR join("\t",$above_below_flag, $chromosome, $start_region, $end_region, $norm_difference, $occup1,$occup2 ), "\n";
		}
    }
	
}



my $size = keys %occupancy;

my $between=$size-$above_counter-$below_counter;
print STDERR "done!\n======================\n";
print STDERR "Upper threshold: ",$threshold1, "\n";
print STDERR "Lower threshold: ",$threshold2, "\n";
print STDERR "======================\n",
	     "$between from $size entries are between 2 thresholds\n",
	     "$above_counter from $size entries saved to $output1\n",
	     "$below_counter from $size entries saved to $output2\n======================\n";
close($OUT1_FHs) or die $!;
close($OUT2_FHs) or die $!;
exit;


#--------------------------------------------------------------------------
# read compressed occupancy file
sub ReadFile {
    my ($in_file, $filename1, $filename2, $col_coords, $col_occup, $col_stdev,$col_relerr, $occupancy_hashref) = @_;
    my $filesize = -s $in_file; #determine file size in bytes
    my $size_counter_step=int($filesize/100);
    $filesize = int($filesize/1048576); # filesize in megabytes

    print STDERR "Reading $in_file file of $filesize MBs. Please wait...\n";

    #read file with by 4kb chanks
    #@coord_occ_array=();
    my $BUFFER_SIZE = 1024*4;
    
	# open compressed occupancy file
	my $inFH;
	if ( $in_file =~ (/.*\.gz$/) ) {
		$inFH = IO::Uncompress::Gunzip->new( $in_file )
		or die "IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip::GunzipError\n";
		
		use constant MIN_COMPRESS_FACTOR => 250;
		my $outer_bytes = -s $in_file;
		my $inner_bytes = get_isize( $in_file );
		$inner_bytes += 4294967296 if( $inner_bytes < $outer_bytes * MIN_COMPRESS_FACTOR );
		$filesize = int($inner_bytes/1048576);
	}
	else { open( $inFH, "<", $in_file ) or die "error: $in_file cannot be opened:$!"; }

    my $buffer = "";
    my $sz_buffer = 0;
    my $timer2 = time();
    # counter for the markers we see
    my $counter = 0;
    
    my $regex_split_newline='\n';
    my $offset=0;
        
    while ((my $n = read($inFH, $buffer, $BUFFER_SIZE)) !=0) {
        if ($n >= $BUFFER_SIZE) {
        $buffer .= <$inFH>;
        }
        my @lines = split(/$regex_split_newline/o, $buffer);
        # process each line in zone file
        foreach my $line (@lines) {
			if ($line =~ /^\D*\t.*/) { next; }
			my @string;
			push (@string, split("\t",$line));
			my $pos = $string[$col_coords]; $pos+=0;
			my $occup = $string[$col_occup]; $occup+=0;
			my $stdev= $string[$col_stdev]; $stdev+=0;
			my $relerr= $string[$col_relerr]; $relerr+=0;

			if ($occup==0) { undef @string; next; }
			
			
			
			$occupancy_hashref->{$pos}->{$filename1} =  $occup;
			$occupancy_hashref->{$pos}->{$filename1}->{'stdev'} = $stdev;
			$occupancy_hashref->{$pos}->{$filename1}->{'relerr'} = $relerr;

			if(! exists ($occupancy_hashref->{$pos}->{$filename2}) ) { $occupancy_hashref->{$pos}->{$filename2}=0; }
			if(! exists ($occupancy_hashref->{$pos}->{$filename2}->{'stdev'}) ) { $occupancy_hashref->{$pos}->{$filename2}->{'stdev'}=0; }
			if(! exists ($occupancy_hashref->{$pos}->{$filename2}->{'relerr'}) ) { $occupancy_hashref->{$pos}->{$filename2}->{'relerr'}=0; }
			undef @string;
			
			if ( $counter % 120000) { }
					else {
						  my $percents = sprintf("%.2f", 100*$offset/(1048576*$filesize));
						  my $processed_fs = sprintf("%.2f", $offset/1048576);
						  my ($progress_seconds,$progress_minutes,$progress_hours,$elapsed_time) = (0,0,0,0);
						  $progress_seconds = time()-$timer2;
						  if ($progress_seconds>=3600) { $progress_hours=floor($progress_seconds/3600); }
						  if ($progress_seconds>=60) { $progress_minutes=floor($progress_seconds/60)-60*$progress_hours; }
						  $progress_seconds=$progress_seconds-3600*$progress_hours-60*$progress_minutes;
						  $elapsed_time=join("", $progress_hours,"h:",$progress_minutes,"m:",$progress_seconds,"s");

						  print STDERR $processed_fs," MBs from $filesize MBs (",$percents,"%) processed in ", $elapsed_time, ".                        \r";
						  }
			$counter++;

        }
        $offset += $n;
        undef @lines;
        $buffer = "";
    }
    
	my ($progress_seconds,$progress_minutes,$progress_hours,$elapsed_time) = (0,0,0,0);
	$progress_seconds = time()-$timer2;
	if ($progress_seconds>=3600) { $progress_hours=floor($progress_seconds/3600); }
	if ($progress_seconds>=60) { $progress_minutes=floor($progress_seconds/60)-60*$progress_hours; }
	$progress_seconds=$progress_seconds-3600*$progress_hours-60*$progress_minutes;
	$elapsed_time=join("", $progress_hours,"h:",$progress_minutes,"m:",$progress_seconds,"s");

    print STDERR "\n",int($offset/1048576), " Mbs processed in ", $elapsed_time, ".\ndone.\n";
	print STDERR "$counter positions added to DB\n";


    close($inFH) or die $!;
    return;
}

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
	if ( ! -f $input2 ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input file $input2: $!\n"
		);
	}
	if ( ! -e $input1 ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input file $input1: $!\n"
		);
	}
	if ( !$output1 )  {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Please specify output regions file name!\n"
		);
	}
	if ( !$output2 )  {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Please specify output regions file name!\n"
		);
	}

}

sub get_isize
{
   my ($file) = @_;

   my $isize_len = 4;

   # create a handle we can seek
   my $FH;
   unless( open( $FH, '<:raw', $file ) )
   {
      die "Failed to open $file: $!";
   }
   my $io;
   my $FD = fileno($FH);
   unless( $io = IO::Handle->new_from_fd( $FD, 'r' ) )
   {
      die "Failed to create new IO::Handle for $FD: $!";
   }

   # seek back from EOF
   unless( $io->IO::Seekable::seek( "-$isize_len", 2 ) ) 
   {
      die "Failed to seek $isize_len from EOF: $!"
   }

   # read from here into mod32_isize
   my $mod32_isize;
   unless( my $bytes_read = $io->read( $mod32_isize, $isize_len ) )
   {
      die "Failed to read $isize_len bytes; read $bytes_read bytes instead: $!";
   }

   # convert mod32 to decimal by unpacking value
   my $dec_isize = unpack( 'V', $mod32_isize );

   return $dec_isize;
}