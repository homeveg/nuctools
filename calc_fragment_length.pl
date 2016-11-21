#!/usr/bin/perl

=head1 NAME

calc_fragment_length.pl -  Calculates frequency of nucleosome-nucleosome distances to determine the nucleosome repeat length

=head1 SYNOPSIS
perl -w calc_fragment_length.pl --input=<in.bed> --output=<filtered.txt> \
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
    --delta | -d                  maximum distance from start of the reference nucleosome to the last in calculations (default: 400)
    --filtering_threshold | -t    remove nucleosome piles above threshold (default: 20)
    --pile | -p                   define minimal pile size (default: 1)
   
   flags:
    --apply_filter | -f           apply --filtering_threshold to the data
    --fix_pile_size | -s          only consider nucleosomes in piles of the defined size (requires -p parameter)

	--help | h                 Help
	
 Example usage:
 
    perl -w calc_fragment_length.pl --input=in.bed.gz --output=filtered.txt.gz --delta=1000 	
	
	OR
	
    perl -w calc_fragment_length.pl -in in.bed.gz -out out.bed.gz -d 100 	
    
=head1 DESCRIPTION
 
=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 calc_fragment_length.pl

 calc_fragment_length.pl -  Estimates mean fragment length for a single-emd sequencing

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
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError) ;

use strict "vars";
use Config;
use Time::localtime;
use Time::Local;
use List::Util qw(first);

# Default parametrs
my $delta = 400;
my $pile = 1;

my $in_file; 
my $out_path1; 
#  Time count Initialisation
my $timer1=time();
my $tm = localtime;
my $start_sec = $tm -> [0];
my $start_min = $tm ->[1];
my $start_hour = $tm ->[2];
my $start_time = time();

my $apply_filter_flag;
my $piles_filtering_threshold=20;
my $fix_pile_size;

# default BED file columns
my $start_col=1;
my $end_col=2;
my $strand_col=5;
my $chromosome_col=0;

my $needsHelp;

#read arguments from command line
my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$in_file,
	'output|out=s'   => \$out_path1,
	
	'delta|d=i' => \$delta,
	'pile|p=i'   => \$pile,
	'filtering_threshold|t=i'   => \$piles_filtering_threshold,
	
	'start_col|sC=s' => \$start_col,
	'end_col|eC=s'   => \$end_col,
	'strand_col|str=s' => \$strand_col,
	'chromosome_col|chr=s'   => \$chromosome_col,

	'fix_pile_size|s' => \$fix_pile_size,
	'apply_filter|f'   => \$apply_filter_flag,
	'help|h'      => \$needsHelp
);

# Check to make sure options are specified correctly and files exist
&check_opts();

# perl -w new_nuc-nuc_distance_filter.pl -input="chr9.bed" -output="nuc-nuc_ch9_filtered.txt" -delta=3000 -filtering_threshold=20 -apply_filter

# Display input parametrs
print STDERR "======================================\n";
print STDERR "Started:\t$start_hour:$start_min:$start_sec\n";
print STDERR "======================================\n";
print STDERR "in file:",$in_file, "\n";
print STDERR "out file:",$out_path1, "\n";
print STDERR "delta: ",$delta, "\n";
print STDERR "pile: $pile\n";
print STDERR "filtering threshold: $piles_filtering_threshold\n";
if ( defined $fix_pile_size) { print STDERR "select only fix pile size: $pile\n"; }
if ( defined $apply_filter_flag) { print STDERR "filter the data: remove all piles above $piles_filtering_threshold\n"; }
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
	or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
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
my $offset=0;
my $not_zero_counter=0;
my $string_counter=0;
my $chr_start;  #first read start
my $chr_end;    # last read end

my (@starts, @ends, @length, %starts_hash, @indexis);
while ((my $n = read($inFH, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <$inFH>;
    }
    my @lines = split(/$regex_split_newline/o, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
	chomp($line);
        my @newline1=split(/\t/, $line);
        my $start_nuc=$newline1[1];
        my $end_nuc=$newline1[2];
        push(@starts, $start_nuc);
        push(@ends, $end_nuc);

        $string_counter++;
	if ($start_nuc>0) {$not_zero_counter++;}
    }
$processed_memory_size += $n;
$offset += $n;
if(int($processed_memory_size/1048576)>= $filesize/10) {
    print STDERR "."; $processed_memory_size=0;
    }
undef @lines;
$buffer = "";
}

close($inFH) or die $!; 
print STDERR $filesize, " Mbs processed in ", time()-$timer2, " seconds.\n$not_zero_counter non zero counts, ",$#starts+1," lines\n\n";

# remove empty strings
$timer2= time();
print STDERR "- cleaning from empty strings (if any)...";
@starts = grep /\S/, @starts;
@ends = grep /\S/, @ends;
print STDERR "done in ", time()-$timer2, " seconds. ",$#starts+1," strings left\n";

# sort nucleosome positions according to a start_nuc
$timer2= time();
print STDERR "- sorting...";
my @sorted_starts = sort {$a <=> $b} @starts;
my @sorted_ends = sort {$a <=> $b} @ends;

print STDERR "done in ", time()-$timer2, " seconds.\n";
$timer2= time();

# remove nucleosomoes without repeat ($pile>1)
my @temp;
if ($pile>1) {
    $timer2= time();
    if($fix_pile_size) {
	print STDERR "- select only pile=$pile...";
    }
    else {
	print STDERR "- removing un-piled nucleosomes...";
    }
    
    my @only_piled;
    my $pile_counter=0;

    for (my $i=1; $i<=$#sorted_starts; $i++) {
	if (!@temp) { push(@temp,$sorted_starts[$i-1]); }
	if ($sorted_starts[$i-1] == $sorted_starts[$i]) {
	    push(@temp,$sorted_starts[$i]);
	    $pile_counter++;
	}
	elsif ($pile_counter < $pile) {
	    undef @temp;
	    $pile_counter=0;
	}
	elsif (($sorted_starts[$i-1] != $sorted_starts[$i]) && ($#temp>0)) {
	    if(($fix_pile_size eq "yes") && ($#temp != $pile)) {
		undef @temp;
	    }
	    else {
		push @only_piled, @temp;
		undef @temp;
	    }
	    $pile_counter=0;
	}

    }
    undef @sorted_starts;
    @sorted_starts = grep /\S/, @only_piled;
    print STDERR "done in ", time()-$timer2, " seconds. ",$#sorted_starts+1," strings left\n";
}

#apply local pile filter: at one position should be less than threshold nucleosom starts
if(! $fix_pile_size ) {
    $timer2= time();
    undef @temp;
    if ($apply_filter_flag) {
	$timer2= time();
	print STDERR "- apply local pile filter: removing reads in the pile above $piles_filtering_threshold ...";
	my @piled_under_threshold;
	my $pile_counter=0;
    
	for (my $i=1; $i<=$#sorted_starts; $i++) {
	    if (!@temp) { push(@temp,$sorted_starts[$i-1]); }
	    if ($sorted_starts[$i-1] == $sorted_starts[$i]) {
		push(@temp,$sorted_starts[$i]);
		$pile_counter++;
	    }
	    elsif ($pile_counter >= $piles_filtering_threshold) {
		push @piled_under_threshold, @temp[0..$piles_filtering_threshold];
		undef @temp;
		$pile_counter=0;
	    }
	    elsif (($sorted_starts[$i-1] != $sorted_starts[$i]) && ($#temp>0)) {
		push @piled_under_threshold, @temp;
		undef @temp;
		$pile_counter=0;
	    }
    
	}
	undef @sorted_starts;
	@sorted_starts = grep /\S/, @piled_under_threshold;
	print STDERR "done in ", time()-$timer2, " seconds. ",$#sorted_starts+1," strings left\n";
    }
}


undef @temp;
if ($pile>1) {
    $timer2= time();
    if($fix_pile_size) { 	print STDERR "- select only pile=$pile...";      }
    else {  	print STDERR "- removing un-piled nucleosomes...";     }
    
    my @only_piled;
    my $pile_counter=0;

    for (my $i=1; $i<=$#sorted_ends; $i++) {
	if (!@temp) { push(@temp,$sorted_ends[$i-1]); }
	if ($sorted_ends[$i-1] == $sorted_ends[$i]) {
	    push(@temp,$sorted_ends[$i]);
	    $pile_counter++;
	}
	elsif ($pile_counter < $pile) {
	    undef @temp;
	    $pile_counter=0;
	}
	elsif (($sorted_ends[$i-1] != $sorted_ends[$i]) && ($#temp>0)) {
	    if(($fix_pile_size eq "yes") && ($#temp != $pile)) {
		undef @temp;
	    }
	    else {
		push @only_piled, @temp;
		undef @temp;
	    }
	    $pile_counter=0;
	}

    }
    undef @sorted_ends;
    @sorted_ends = grep /\S/, @only_piled;
    print STDERR "done in ", time()-$timer2, " seconds. ",$#sorted_ends+1," strings left\n";
}

#apply local pile filter: at one position should be less than threshold nucleosom starts
if(! $fix_pile_size ) {
    $timer2= time();
    undef @temp;
    if ($apply_filter_flag) {
	$timer2= time();
	print STDERR "- apply local pile filter: removing reads in the pile above $piles_filtering_threshold ...";
	my @piled_under_threshold;
	my $pile_counter=0;
    
	for (my $i=1; $i<=$#sorted_ends; $i++) {
	    if (!@temp) { push(@temp,$sorted_ends[$i-1]); }
	    if ($sorted_ends[$i-1] == $sorted_ends[$i]) {
		push(@temp,$sorted_ends[$i]);
		$pile_counter++;
	    }
	    elsif ($pile_counter >= $piles_filtering_threshold) {
		push @piled_under_threshold, @temp[0..$piles_filtering_threshold];
		undef @temp;
		$pile_counter=0;
	    }
	    elsif (($sorted_ends[$i-1] != $sorted_ends[$i]) && ($#temp>0)) {
		push @piled_under_threshold, @temp;
		undef @temp;
		$pile_counter=0;
	    }
    
	}
	undef @sorted_ends;
	@sorted_ends = grep /\S/, @piled_under_threshold;
	print STDERR "done in ", time()-$timer2, " seconds. ",$#sorted_ends+1," strings left\n";
    }
}

print STDERR "Step 2: Start analysis...\n";
$timer2= time();    
my @output_array = (0) x ($delta+1); #initialize array of 0 of $delta+1 size
my $first_itteration=0;
my $counter_step=int($#sorted_starts/100);
my $progress_counter=$counter_step;

for (my $i=0; $i<$#sorted_starts; $i++) {
    #read read start
    my $nuc_start=$sorted_starts[$i];
    # calcualte maximum index shift
    my $max_delta_index=5*$delta;
    # check if incremented index exceeds reads array length and correct it if necessary
    if ($i+$max_delta_index>=$#sorted_starts) { $max_delta_index = $#sorted_starts-$i; }
    for (my $n=1; $n<=$max_delta_index ; $n++) {
	if (!$sorted_starts[$i+$n]) { last; }
        my $nuc_plus_end=$sorted_ends[$i+$n];
        my $delta_nuc_starts = $nuc_plus_end-$nuc_start;
	if (($delta_nuc_starts>$delta) or ($delta_nuc_starts<0)) { last; }
        else { $output_array[$delta_nuc_starts]++; }
    }
    #increment counter to display work progress
    if ( $progress_counter == $i ) {
	if ($first_itteration==0) { my $approx_duration = (time()-$timer2)*100/60; print STDERR "work will be finished in approximately $approx_duration minutes\n"; $first_itteration=1; }
	print STDERR ".";
	#print STDERR join (" ", @output_array[0..5],"..",@output_array[145..150]), "\n";
	$progress_counter+=$counter_step;
    }
}
print STDERR "\ndone\n"; 

@output_array = grep /\S/, @output_array;

print STDERR "- saving resuts to $out_path1...";

# open pipe to Gzip or open text file for writing
my $out_file = $out_path1;
$out_file =~ s/(.*)\.gz$/$1/;
my $gz_out_file = $out_file.".gz";
my $OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
print $OUT_FHs join("\n", @output_array);
close (OCCUP_File);
print STDERR "done\n\n";
my $fragment_length_estimate = first { $output_array[$_] >0 } 0..$#output_array;
print STDERR "Estimated fragment length is $fragment_length_estimate \n\n";

my $out_path_fl;
if ($out_path1 =~/\//) {
    $out_path_fl=$out_path1;
    $out_path_fl =~s/(.*\/)(.*)$/$1Fragment_length.$2/;
}
else { $out_path_fl= "Fragment_length.".$out_path1; }

print STDERR "saving fragment length to to a file $out_path_fl...";
open (OCCUP_File, ">$out_path_fl") or die "can't open file $out_path_fl for writting: $!";
print OCCUP_File $fragment_length_estimate,"\n";
close (OCCUP_File);
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

undef @output_array;
exit;


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

#===========================================================#
sub median {
	my $median;
	my @data = sort {$a <=> $b } @_;
	if (even_odd($#data+1) eq "ODD") {
	$median = $data[$#data/2];
	}
	else {
	$median = Average($data[$#data/2],$data[($#data/2)+1]);
	}
	return($median);
}
#===========================================================#
sub even_odd {
	if (int($_[0]/2) == $_[0]/2) { return "EVEN"; }
	else {return "ODD";}
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
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input BED file $in_file: $!\n"
		);
	}
	if (!$out_path1 ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify output file name\n"
		);
	}

}
