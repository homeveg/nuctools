#!/usr/bin/perl

=head1 NAME

nucleosome_repeat_length.pl -  Calculates frequency of nucleosome-nucleosome distances to determine the nucleosome repeat length

=head1 SYNOPSIS

perl -w nucleosome_repeat_length.pl --input=<in.bed> --output=<filtered.txt> \
 [--delta=<N> --apply_filter --filtering_threshold=<N> --pile=<N> --fix_pile_size ] \
 [--chromosome_col=<column Nr.> --start_col=<column Nr.> --end_col=<column Nr.> --strand_col=<column Nr.> --help]


 Required arguments:
    --input | -in      path to input BED or BED.GZ file
    --output | -out    output table file name
	
 Options:
 
  define column numbers in the input BED file (Nr. of the very first column is 0):
    --start_col | -sC            read start column Nr. (default: -s 1)
    --end_col | -eC              read end column Nr. (default: -e 2)
    --strand_col | -str          strand column Nr. (default: -str 5)
    --chromosome_col | -chrC     chromosome column Nr. (default: -chr 0)

   parameters with default values:
    --delta | -d                  maximum distance from start of the reference nucleosome to the last in calculations (default: 1500)
    --pile_delta | -pD            maximum distance between adjacent nucleosome starts to consider as one pile (default: 5)
    --filtering_threshold | -t    remove nucleosome piles above threshold (default: 20)
    --pile | -p                   define minimal pile size (default: 1)
	--MaxNr | -m                  set maximum number of adjacent reads to analyze (default: 1000000 )

   flags:
    --apply_filter | -f           apply --filtering_threshold to the data
	--detectPiles | -dP           auto-detect piles of considerable size and apply --filtering_threshold to the data
    --fix_pile_size | -s          only consider nucleosomes in piles of the defined size (requires -p parameter)
	
   alignment settings (default behavior: align mid point):
    --useStrand | -uS             use strand information (for single-end sequencing reads)
	--alignStarts | -aS           use read start to align nucleosomes. Please note: --alignStarts, --alignEnds and --useStrand flags are mutually exclusive!
	--alignEnds | -aE             use read ends to align nucleosomes. Please note: --alignStarts, --alignEnds and --useStrand flags are mutually exclusive!

	--help | h                    Help
	
 Example usage:
 
    perl -w nucleosome_repeat_length.pl --input=in.bed.gz --output=filtered.txt.gz --delta=1000 	
	
	OR
	
    perl -w nucleosome_repeat_length.pl -in in.bed.gz -out out.bed.gz -d 100 	
    
=head1 DESCRIPTION
 
=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 nucleosome_repeat_length.pl

 nucleosome_repeat_length.pl -  Calculates frequency of nucleosome-nucleosome distances to determine the nucleosome repeat length

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
use POSIX;
use File::Basename;
use Time::localtime;
use Time::Local;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }




# Default parameters
my $delta = 1500;
my $pile = 1;
my $piles_filtering_threshold=20;
my $pile_delta = 5;
my $detectPiles;
my $in_file; 
my $out_path1;

# default BED file columns
my $start_col=1;
my $end_col=2;
my $strand_col=3;
my $chromosome_col=0;

#  Time count Initialization
my $timer1=time();
my $tm = localtime;
my $start_sec = $tm -> [0];
my $start_min = $tm ->[1];
my $start_hour = $tm ->[2];
my $start_time = time();

# initialize flags
my $apply_filter_flag;
my $fix_pile_size;
my $MaxNr=1000000;

my $needsHelp;
my $useGZ;
my $useStrand;
my $alignStarts;
my $alignEnds;

#read arguments from command line
my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$in_file,
	'output|out=s'   => \$out_path1,
	
	'delta|d=i' => \$delta,
	'pile|p=i'   => \$pile,
	'pile_delta|pD=i'   => \$pile_delta,
	'filtering_threshold|t=i'   => \$piles_filtering_threshold,
	'MaxNr|m=i' => \$MaxNr,

	'start_col|sC=s' => \$start_col,
	'end_col|eC=s' => \$end_col,
	'strand_col|str=s' => \$strand_col,
	'chromosome_col|chr=s'   => \$chromosome_col,

	'fix_pile_size|s' => \$fix_pile_size,
	'apply_filter|f'  => \$apply_filter_flag,
	'detectPiles|dP'  => \$detectPiles,
	'useStrand|uS' => \$useStrand,
	
	'alignStarts|aS' => \$alignStarts,
	'alignEnds|aE' => \$alignEnds,

	'gzip|z' => \$useGZ,

	'help|h'      => \$needsHelp
);

# Check to make sure options are specified correctly and files exist
&check_opts();

# check if GZIP is loaded
if ((!$ModuleGzipIsLoaded) or (!$ModuleGunzipIsLoaded))  {
	print STDERR "Can't work with GZIP: IO::Compress::Gzip is not on PATH\n";
	exit;
}
elsif (($ModuleGzipIsLoaded) and ($ModuleGunzipIsLoaded) ) {
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

# Display input parameters
print STDERR "======================================\n";
print STDERR "Started:\t$start_hour:$start_min:$start_sec\n";
print STDERR "======================================\n";
print STDERR "in file:",$in_file, "\n";
print STDERR "out file:",$out_path1, "\n";
print STDERR "delta: ",$delta, "\n";
print STDERR "pile: $pile\n";
print STDERR "filtering threshold: $piles_filtering_threshold\n";
if ( defined $fix_pile_size) { print STDERR "select only fix pile size: $pile\n"; }
print STDERR "pile delta: ",$pile_delta, "\n";
if ( defined $apply_filter_flag) { print STDERR "filter the data: remove all piles above $piles_filtering_threshold\n"; }
if ( defined $detectPiles) { print STDERR "piles auto-detection has been activated. Piles filtering threshold set to $piles_filtering_threshold\n";}
if ( defined $alignStarts ) { print STDERR "single-end analysis mode: use read start to align nucleosomes\n";
							  print STDERR "Please note: --alignStarts, --alignEnds and --useStrand flags are mutually exclusive!\n";
							  print STDERR "disabling --alignEnds and --useStrand flags\n";
							  undef $alignEnds; undef $useStrand;}
if ( defined $alignEnds ) { print STDERR "single-end analysis mode: use read ends to align nucleosomes\n";
							  print STDERR "Please note: --alignStarts, --alignEnds and --useStrand flags are mutually exclusive!\n";
							  print STDERR "disabling --alignStarts and --useStrand flags\n";
							  undef $alignStarts; undef $useStrand;}
if ( defined $useStrand ) { print STDERR "analyzing single-end reads: use strand from column $strand_col\n"; }

else { print STDERR "analyzing paired-end reads: ignore strand\n"; }
print STDERR "Nr of reads to limit fragment length calculation: $MaxNr\n";
print STDERR "======================================\n";


my @occ_array=();

#read first file with occupanicies
#read file with by 4kb chanks

print STDERR "Step 1 (",time()-$timer1," sec.): reading $in_file file...\n";

@occ_array=();
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

my $regex_split_tab='\t';
my $regex_split_newline='\n';

my $filesize = -s $in_file; #determine file size in bytes
my $size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "- reading nucleosome start position column from $in_file ($filesize MBs). Please wait...\n";

my $processed_memory_size = 0;
my $total_mbytes_loaded=0;
my $offset=0;
my $not_zero_counter=0;
my $string_counter=0;
my $chr_start;  #first read start
my $chr_end;    # last read end
my $cancel_load;

my (%hash);
while ((my $n = read($inFH, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <$inFH>;
    }
    my @lines = split(/$regex_split_newline/o, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
		chomp($line);
		my @row=split(/\t/, $line);
		my $start_nuc=$row[$start_col];
		my $end_nuc=$row[$end_col];

		my $strand;
		if ($useStrand) { $strand = $row[$strand_col] eq '+' ? 'plus' : 'minus' ; }
		else { $strand = 'plus'; }
		
		if ($useStrand) {
			$hash{$string_counter}{$strand}{start}=$start_nuc;
			$hash{$string_counter}{$strand}{end}=$end_nuc;
		}
		elsif ($alignStarts) {
			$hash{$string_counter}{start}=$start_nuc;
		}
		elsif ($alignEnds) {
			$hash{$string_counter}{end}=$start_nuc;
		}
		$hash{$string_counter}{mid}=min($end_nuc,$start_nuc)+int(abs($end_nuc-$start_nuc)/2);

		$string_counter++;
		if ($start_nuc>0) {$not_zero_counter++;}
		if ( $string_counter == $MaxNr ) {
			print STDERR "reach read number limit $MaxNr. Proceeding to the next steps...\n";
			$cancel_load="yes";
			last;
		}
    }
	$processed_memory_size += $n;
	$total_mbytes_loaded += $n;
	$offset += $n;
	if(int($processed_memory_size/1048576)>= $filesize/10) {
		print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.            \r";
		$processed_memory_size=0;
    }
	undef @lines;
	$buffer = "";
	if($cancel_load) {last;}

}

close($inFH) or die $!;
$total_mbytes_loaded = sprintf "%.2f", $processed_memory_size/1048576;

print STDERR $total_mbytes_loaded, " Mbs processed in ", time()-$timer2, " seconds.\n$not_zero_counter non zero counts\n\n";
if ( $string_counter < $MaxNr ) {
	print STDERR "file $in_file loaded to 100%\n\n"
}
else
{
	my $rounded = sprintf "%.2f", 100*$total_mbytes_loaded/$filesize;
	print STDERR "file $in_file loaded to ", $rounded  , "%\n\n";
}
# sort nucleosome positions according to a start_nuc
$timer2= time();
print STDERR "sorting...";

my $data = \%hash;
my @sorted_array;

my @sorted_hash = sort {
	$data->{$a}{mid} <=> $data->{$b}{mid}
	or
	$data->{$a}{mid} cmp $data->{$b}{mid}
} keys %$data;
	
print STDERR "done in ", time()-$timer2, " seconds.\n\n";

my (@sorted_plus_starts, @sorted_minus_starts, @sorted_mids, @sorted_starts, @sorted_ends );
foreach my $k (@sorted_hash) {
    #print "$k: $data->{$k}{mid}\n";
	if ($useStrand) {
		if($data->{$k}{plus}{start}) { push (@sorted_plus_starts, $data->{$k}{plus}{start}); }
		if($data->{$k}{minus}{start}) { push (@sorted_minus_starts, $data->{$k}{minus}{start}); }
	}
	elsif ($alignStarts) { push (@sorted_starts, $data->{$k}{start}); }
	elsif ($alignEnds) { push (@sorted_ends, $data->{$k}{end}); }
	else { if($data->{$k}{mid}) { push (@sorted_mids, $data->{$k}{mid}); }
	}
}


if ($useStrand) {
	print STDERR join("\t", "plus starts: $#sorted_plus_starts",  "minus starts: $#sorted_minus_starts"),"\n";
	for (my $i=0; $i<=10; $i++) {
		print STDERR join("\t", $sorted_plus_starts[$i], $sorted_minus_starts[$i]),"\n";
	}
} elsif ($alignStarts) {
	print STDERR "starts: $#sorted_starts\n";
	for (my $i=0; $i<=10; $i++) {
		print STDERR "$sorted_starts[$i]\n";
	}
} elsif ($alignEnds) {
	print STDERR "ends: $#sorted_ends\n";
	for (my $i=0; $i<=10; $i++) {
		print STDERR "$sorted_ends[$i]\n";
	}
} else {
	print STDERR "mids: $#sorted_mids\n";
	for (my $i=0; $i<=10; $i++) {
		print STDERR "$sorted_mids[$i]\n";
	}
}

print STDERR "done in ", time()-$timer2, " seconds.\n\n";
$timer2= time();

# check for piles of considerable size
my @piles;
if ($useStrand) {
	@piles=check_piles($piles_filtering_threshold, floor($pile_delta/2), @sorted_plus_starts);
} elsif ($alignStarts) {
	@piles=check_piles($piles_filtering_threshold, floor($pile_delta/2), @sorted_starts);
} elsif ($alignEnds) {
	@piles=check_piles($piles_filtering_threshold, floor($pile_delta/2), @sorted_ends);
} else{
	@piles=check_piles($piles_filtering_threshold, floor($pile_delta/2), @sorted_mids);
}

my $max_pile=max(@piles);
print STDERR $#piles+1," piles above $piles_filtering_threshold has been detected.\n",
"Max pile size: $max_pile\n",
"Average pile size: ",floor(average(\@piles))," \x{B1} ", floor( stdev(\@piles) ), "\n",
"Median pile size: ",median(@piles),"\n";

if( ($max_pile>=$piles_filtering_threshold) && ( floor( stdev(\@piles) ) > floor(average(\@piles)) ) ){
	print STDERR "\n================================================================\n",
	"        WARNING: low complexity regions has been detected!\n";
	
	if($detectPiles) {
		$apply_filter_flag=1;
		$piles_filtering_threshold=$piles_filtering_threshold;
		print STDERR "        automatically set piles filtering threshold to ",$piles_filtering_threshold,"\n";
	}
	else {
		print STDERR "        consider setting up piles filtering threshold \nor use piles auto-detection (flag: --detectPiles )\n";
	}
	print STDERR "================================================================\n\n";
}

# remove nucleosomoes without repeats ($pile>1)
if ($pile>1) {
	print STDERR "remove nucleosomoes without repeats\n";
	my @temp;

	if ($useStrand) {
		@temp = remove_unpiled($pile, $fix_pile_size, floor($pile_delta/2), @sorted_minus_starts);
		undef @sorted_minus_starts;
		@sorted_minus_starts = @temp;
		undef @temp;
		
		@temp = remove_unpiled($pile, $fix_pile_size, floor($pile_delta/2), @sorted_plus_starts);
		undef @sorted_plus_starts;
		@sorted_plus_starts = @temp;
		undef @temp;
	}
	else {
		@temp = remove_unpiled($pile, $fix_pile_size, floor($pile_delta/2), @sorted_mids);
		undef @sorted_mids;
		@sorted_mids = @temp;
		undef @temp;
	}

}

if ($apply_filter_flag) {
	print STDERR "remove piles above $piles_filtering_threshold\n";
	my @temp;
		
	if ($useStrand) {
		@temp = filter_by_threshold($piles_filtering_threshold, floor($pile_delta/2), @sorted_plus_starts);
		undef @sorted_plus_starts;
		@sorted_plus_starts = @temp;
		undef @temp;
		
		@temp = filter_by_threshold($piles_filtering_threshold, floor($pile_delta/2), @sorted_minus_starts);
		undef @sorted_minus_starts;
		@sorted_minus_starts = @temp;
		undef @temp;
	} elsif ($alignStarts) {
		@temp = filter_by_threshold($piles_filtering_threshold, floor($pile_delta/2), @sorted_starts);
		undef @sorted_starts;
		@sorted_starts = @temp;
		undef @temp;		
	} elsif ($alignEnds) {
		@temp = filter_by_threshold($piles_filtering_threshold, floor($pile_delta/2), @sorted_ends);
		undef @sorted_ends;
		@sorted_ends = @temp;
		undef @temp;		
	} else {
		@temp = filter_by_threshold($piles_filtering_threshold, floor($pile_delta/2), @sorted_mids);
		undef @sorted_mids;
		@sorted_mids = @temp;
		undef @temp;	
	}

}

if ($fix_pile_size ) {
	print STDERR "select only piles of size $pile\n";
	my @temp;
	if ($useStrand) {
		my @temp = local_pile_filter($pile, floor($pile_delta/2), @sorted_plus_starts);
		undef @sorted_plus_starts;
		@sorted_plus_starts = @temp;
		undef @temp;
		
		@temp = local_pile_filter($pile, floor($pile_delta/2), @sorted_minus_starts);
		undef @sorted_minus_starts;
		@sorted_minus_starts = @temp;
		undef @temp;
	} elsif ($alignStarts) {
		@temp = local_pile_filter($pile, floor($pile_delta/2), @sorted_starts);
		undef @sorted_starts;
		@sorted_starts = @temp;
		undef @temp;		
	} elsif ($alignEnds) {
		@temp = local_pile_filter($pile, floor($pile_delta/2), @sorted_ends);
		undef @sorted_ends;
		@sorted_ends = @temp;
		undef @temp;
	} else {
		@temp = local_pile_filter($pile, floor($pile_delta/2), @sorted_mids);
		undef @sorted_mids;
		@sorted_mids = @temp;
		undef @temp;
	}
}

print STDERR "\nCalculating phasograms...\n";
$timer2= time();
my (@output_plus_array,@output_minus_array,@output_mids_array,@output_starts_array,@output_ends_array);
if ($useStrand) {
	@output_minus_array=phasogram(\@sorted_minus_starts, $delta);
	@output_plus_array=phasogram(\@sorted_plus_starts, $delta);
	} elsif ($alignStarts) {
	@output_starts_array=phasogram(\@sorted_starts, $delta);
	} elsif ($alignEnds) {
	@output_ends_array=phasogram(\@sorted_ends, $delta);
	} else {
	@output_mids_array=phasogram(\@sorted_mids, $delta);
}

print STDERR "saving results to $out_path1...";
# open pipe to text file for writing
open my $OUT_FHs, '>', $out_path1 or die "Can't open $out_path1 for writing; $!\n";
my $maxind;
if ($useStrand) { $maxind=$#output_plus_array;}
elsif ($alignStarts) { $maxind=$#output_starts_array;}
elsif ($alignEnds) { $maxind=$#output_ends_array;}
else { 	$maxind=$#output_mids_array; }

for (my $ind=0; $ind<=$maxind; $ind++) {
	if ($useStrand) { print $OUT_FHs join("\t", $output_plus_array[$ind], $output_minus_array[$ind]),"\n"; }
    elsif ($alignStarts) { print $OUT_FHs $output_starts_array[$ind],"\n"; }
	elsif ($alignEnds) { print $OUT_FHs $output_ends_array[$ind],"\n"; }
	else { print $OUT_FHs $output_mids_array[$ind],"\n"; }
}
close ($OUT_FHs);
print STDERR "done\n";

$tm = localtime;
my $stop_sec = $tm -> [0];
my $stop_min = $tm ->[1];
my $stop_hour = $tm ->[2];
my $message;

my $duration = time()-$start_time;
print STDERR "======================================\n";
$message = "\nFinished:\t$stop_hour:$stop_min:$stop_sec\nduration:\t$duration sec.\n\n";
print STDERR "$message\nBye!\n";

undef @output_plus_array;
undef @output_minus_array;
undef @output_mids_array;
undef @output_starts_array;
undef @output_ends_array;

exit;


#################################################################################

sub phasogram {
	my ($sorted_starts_ref, $delta) = @_;
	
	my @output_array = (0) x ($delta+1); #initialize array of 0 of $delta+1 size
	my $first_itteration=0;
	my @sorted_starts = @{ $sorted_starts_ref };

	my $counter_step=int($#sorted_starts/100);
	my $progress_counter=$counter_step;

	print STDERR "calculating distances between adjacent nucleosomes (in the region from up to $delta bases away from the origin)...\n";
	
	#print STDERR join("\t", "SIndx", "EIndx", "summ", "nuc_start", "nuc_end", "delta", "total starts", "total ends"), "\n";

	for (my $i=0; $i<$#sorted_starts; $i++) {
		#read read start
		my $nuc_start=$sorted_starts[$i];
		# calcualte maximum index shift
		my $max_delta_index=5*$delta;
		# check if incremented index exceeds reads array length and correct it if necessary
		if ($i+$max_delta_index>=$#sorted_starts) { $max_delta_index = $#sorted_starts-$i; }
		my @other_starts=@sorted_starts[$i..($i+$max_delta_index)];

		for (my $n=0; $n<=$max_delta_index ; $n++) {
			if (!$sorted_starts[$i+$n]) { last; }
			# remove ends 
			for (my $j=0; $j<=$#other_starts; $j++) { 
				if ( $other_starts[0] > $nuc_start ) { last; }
				else { shift @other_starts; }
			}
			my $other_nuc_start=$other_starts[$n];
			my $delta_nuc_starts = abs($other_nuc_start-$nuc_start);
			if ($delta_nuc_starts>$delta) {
				#print STDERR join("\t", "- $i", $n, $i+$n, $nuc_start, $other_nuc_start, $delta_nuc_starts, $#sorted_starts, $#other_starts), "\n";
				$output_array[$delta_nuc_starts]++;
				last;
			}
			else {
				$output_array[$delta_nuc_starts]++;
				#print STDERR join("\t", $i, $n, $i+$n, $nuc_start, $other_nuc_start, $delta_nuc_starts, $#sorted_starts, $#other_starts), "\n";
			}
		}
		undef @other_starts;
		#increment counter to display work progress
		if ( $progress_counter == $i ) {
			$progress_counter+=$counter_step;
			# last; # enable for testing
            my $rounded = sprintf "%.2f", 100*$i/$#sorted_starts;
            print STDERR $rounded,"% done             \r";
		}
	}
	my @results = grep /\S/, @output_array;
	print STDERR "\ndone in ", time()-$timer2, " seconds. ",$#results+1," strings left\n\n";
	return(@results);
}

#--------------- clean line endings ----------------------------
sub clean {

my $text = shift;

$text =~ s/\r//g;
$text =~ s/\n//g;
return $text;
}

#-------------- determine maxsimum value from the array ---------
sub max {
  my $max = $_[0];
  for ( @_[ 1..$#_ ] ) { $max = $_ if $_ > $max; }
  return($max);
}

#-------------- determine maxsimum value from the array ---------
sub min {
  my $min = $_[0];
  for ( @_[ 1..$#_ ] ) { $min = $_ if $_ < $min; }
  return($min);
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array!\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
sub median {
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) {
        return $vals[int($len/2)];
    }
    else {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
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
	if ( ! -e $in_file ) {
		if ( ! -e "$in_file.gz" ) {
			pod2usage(
				-exitval => 2,
				-verbose => 1,
				-message => "Cannot find input BED file $in_file or $in_file.gz: $!\n"
			);
		}
	}
	if (!$out_path1 ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify output file name\n"
		);
	}

}
#============================================================================
sub remove_unpiled {
	my ($pile, $fix_pile_size, $pile_delta, @sorted_coords) = @_;
	my @only_piled = ();
	my @temp;

	$timer2= time();
	if($fix_pile_size) {
	print STDERR "- select only pile=$pile...";
	}
	else {
	print STDERR "- removing un-piled nucleosomes...";
	}
	
	my $pile_counter=1;

	for (my $i=1; $i<=$#sorted_coords; $i++) {
		if (!@temp) { push(@temp,$sorted_coords[$i-1]); }
		if ( ($sorted_coords[$i] >= $sorted_coords[$i-1] ) && ($sorted_coords[$i] <= $sorted_coords[$i-1] + $pile_delta) ) {
			push(@temp,$sorted_coords[$i]);
			$pile_counter++;
		}
		elsif ($pile_counter < $pile) {
			undef @temp;
			$pile_counter=0;
		}
		elsif ( not ( ($sorted_coords[$i] >= $sorted_coords[$i-1] ) && ($sorted_coords[$i] <= $sorted_coords[$i-1] + $pile_delta)) && ($#temp>0) ) {
			if(($fix_pile_size) && ($#temp != $pile)) { undef @temp; }
			else {
				push @only_piled, @temp;
				undef @temp;
			}
			$pile_counter=0;
		}

	}
	my @results = grep /\S/, @only_piled;
	print STDERR "done in ", time()-$timer2, " seconds. ",$#results+1," strings left\n\n";
	return(@results);
	
}

#============================================================================
sub filter_by_threshold {
	my ($piles_filtering_threshold, @sorted_coords) = @_;
	print STDERR "- apply local pile filter: removing reads in the pile above $piles_filtering_threshold ...";
	my @piled_under_threshold;
	my $pile_counter=0;
	my @temp;
	
	for (my $i=1; $i<=$#sorted_coords; $i++) {
		if (!@temp) { push(@temp,$sorted_coords[$i-1]); }
		if ( ($sorted_coords[$i] >= $sorted_coords[$i-1] ) && ($sorted_coords[$i] <= $sorted_coords[$i-1] + $pile_delta) ) {
			push(@temp,$sorted_coords[$i]);
			$pile_counter++;
		} elsif ($pile_counter >= $piles_filtering_threshold) {
			push @piled_under_threshold, @temp[0..$piles_filtering_threshold];
			undef @temp;
			$pile_counter=0;
		} elsif ($#temp>0) {
			push @piled_under_threshold, @temp;
			undef @temp;
			$pile_counter=0;
		}
	
	}
	my @results = grep /\S/, @piled_under_threshold;
	print STDERR "done in ", time()-$timer2, " seconds. ",$#results+1," strings left\n\n";
	return(@results);
}

#============================================================================
sub check_piles {
	my ($piles_filtering_threshold, $pile_delta, @sorted_coords) = @_;
	print STDERR "- check piles above $piles_filtering_threshold ...";
	my $pile_counter=0;
	my @piles;
	my @temp;
	
	for (my $i=1; $i<=$#sorted_coords; $i++) {
		if (!@temp) { push(@temp,$sorted_coords[$i-1]); }
		if ( ($sorted_coords[$i] >= $sorted_coords[$i-1]-$pile_delta ) && ($sorted_coords[$i] <= $sorted_coords[$i-1]+$pile_delta) ) {
			push(@temp,$sorted_coords[$i]);
			$pile_counter++;
		} elsif ($pile_counter >= $piles_filtering_threshold) {
			push @piles, $pile_counter;
			
			undef @temp;
			$pile_counter=0;
		} 
	
	}
	print STDERR "done in ", time()-$timer2, " seconds. \n\n";
	return(@piles);
}

#============================================================================
sub local_pile_filter {
	my ($piles_filtering_threshold, @sorted_coords) = @_;
	my @only_piled = ();
	my @temp;

	print STDERR "- apply local pile filter: removing reads in the pile above $piles_filtering_threshold ...";
	my @piled_under_threshold;
	my $pile_counter=0;
	
	for (my $i=1; $i<=$#sorted_coords; $i++) {
		if (!@temp) { push(@temp,$sorted_coords[$i-1]); }
		if ( ($sorted_coords[$i] >= $sorted_coords[$i-1] ) && ($sorted_coords[$i] <= $sorted_coords[$i-1] + $pile_delta) ) {
			push(@temp,$sorted_coords[$i]);
			$pile_counter++;
		} elsif ($pile_counter >= $piles_filtering_threshold) {
			push @piled_under_threshold, @temp[0..$piles_filtering_threshold];
			undef @temp;
			$pile_counter=0;
		} elsif ($#temp>0) {
			push @piled_under_threshold, @temp;
			undef @temp;
			$pile_counter=0;
		}
	
	}
	my @results = grep /\S/, @piled_under_threshold;
	print STDERR "done in ", time()-$timer2, " seconds. ",$#results+1," strings left\n\n";
	return(@results);

}

#============================================================================
sub sum_arrays_by_row {
	my @arrays = @_;
	
	my $length = $#{ $arrays[0] };
	my @out;
	for my $i (0 .. $length) {
	  my $accumulator = 0;
	  for my $array (@arrays) {
	    $accumulator += $array->[$i];
	  }
	  push @out, $accumulator;
	}
	return(@out);
}