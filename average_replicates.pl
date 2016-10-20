#!/usr/bin/perl

=head1 NAME

average_replicates.pl - Calculates the average occupancy profile based on several replicate occupancy profiles 

=head1 SYNOPSIS

perl -w average_replicates.pl --input=<path to working dir> --output=<path to results file> --coordsCol=0 --occupCol=1 --pattern="occ.gz" --printData [--help] 

 Required arguments:
    --input | -in      path to directory with aggregate profiles
    --output | -out    output table file name
	
 Options:
   define column numbers in the input occupancy files (Nr. of the very first column is 0):
    --coordsCol | -cC  chromosome coordinates column Nr. (default: -cC 1)
    --occupCol | -oC   cumulative occupancy column Nr. (default: -oC 0)

   additional parameters
	--pattern | -p     occupancy profile file name extension template (default: occ.gz)
    --printData | -d   print all input occupancy columns to the output file
    --help | -h        Help
    
 Example usage:
    perl -w average_replicates.pl --input=/mnt/hdd01/myWD --output=/mnt/hdd01/myWD/occup_tab.txt --pattern="occ.gz" --coordsCol=0 --occupCol=1 --printData
	
	OR
	
    perl -w average_replicates.pl -in /mnt/hdd01/myWD -out /mnt/hdd01/myWD/occup_tab.txt -p "occ.gz" -cC 0 -oC 1 -d
    
=head1 DESCRIPTION
 
=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 average_replicates.pl

 extract_chr_bed.pl calculates the average occupancy profile and standard deviation based on several replicate occupancy profiles from the working directory and save resulting table, including input occupancy data for individual files. Input *.occ files can be flat or compressed. Resulting extended occupancy file will be saved compressed 

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

use Time::localtime;
use Time::Local;
use File::Basename;
use List::Util qw(sum);

use strict 'vars';
use Getopt::Long;
use Pod::Usage;
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError) ;


my (%occupancy,%NormFactors);

my $wd;
my $output;
my $filename_pattern='occ.gz';

my $coordsCol=0;
my $occupCol=1;
my $addData;

my $needsHelp;

my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$wd,
	'output|out=s'   => \$output,
	'pattern|p=s' => \$filename_pattern,

	'coordsCol|cC=s' => \$coordsCol,
	'occupCol|oC=s' => \$occupCol,
	
	'printData|d' => \$addData,
	'help|h'      => \$needsHelp
);

# set flags
$addData = $addData ? "yes" : "no";
# Check to make sure options are specified correctly and files exist
&check_opts();


my $tm = localtime;
print STDERR "-----------------------\n",
join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",
join(":",$tm -> [2],$tm -> [1],$tm -> [0]),
"\n-----------------------\n";

#check if folder exists

opendir(DIR, "$wd") or die $!;
my @all_files = readdir(DIR);
closedir(DIR);
my (@names,@files);

foreach my $file (sort @all_files){
  if ($file =~ m/.*\.$filename_pattern$/){
	push(@files, $file);
	my $filename = basename($file,  "\.$filename_pattern");
	push(@names, $filename);
	print STDERR "process $filename...\n";
	$NormFactors{$filename} = ReadFile("$wd/$file", $filename, $coordsCol, $occupCol, \%occupancy);
	}
}

print STDERR "calculating StDev, Variance and average.\nResults will be saved to $output\n";

# open pipe to Gzip or open text file for writing
my $out_file = $output;
$out_file =~ s/(.*)\.gz$/$1/;
my $gz_out_file = $out_file.".gz";
my $OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";

my $size = keys %occupancy;
print $OUT_FHs join("\t","position","Mean","stdev","Rel.Error");
my $header=0;
my $total_counts = keys %occupancy;
print STDERR "processing $total_counts entries...\n";
my $work_progress_step = int($total_counts/10);
my $current_progress = $work_progress_step;
my $j=0;
for my $position ( sort {$a<=>$b} keys %occupancy) {
    if($current_progress == $j) {print STDERR "$current_progress from $total_counts...\n";$current_progress+=$work_progress_step;}
    $j++;
    my @temp; my $coord;
    for my $name ( sort keys %{ $occupancy{$position} }) {
	if ( ($header==0) && ($addData eq "yes") ){ print $OUT_FHs "\t$name"; }
	my $occup;
	if ($NormFactors{$name}==0) {
	    print STDERR "index:$j\tname:$name\tposition:$position\tnorm factor:$NormFactors{$name}\toccupancy:$occupancy{$position}{$name}\tratio:NA\n";
	}
	else {
	    $occup=$occupancy{$position}{$name}/$NormFactors{$name};

	}

	push(@temp, $occup);
    }
    if ($header == 0) { print $OUT_FHs "\n"; }
    $header=1;
    my $Mean=calcMean(@temp);
    my ($stdev,$variance)=calcStdDev(@temp);
    my $rel_err;
    if ($Mean!=0) {  $rel_err=$stdev/$Mean; }
    else {$rel_err=0;}

    if ($addData eq "yes") { print $OUT_FHs join("\t",$position,$Mean,$stdev,$rel_err,@temp),"\n";   }
    else {print $OUT_FHs join("\t",$position,$Mean,$stdev,$rel_err),"\n";   }
    undef @temp;
}
print STDERR "done!\n";
$tm = localtime;
print STDERR "-----------------------\n",join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),"\n-----------------------\n";
close($OUT_FHs) or die $!;
exit();



###################################


#------------------ Read specified column --------------------
sub Read_column {
    my ($column_number, $array_ref) = @_;
    my @array = @{$array_ref};

    my (@column, @string);
    # read column of interest to the memory
    for(my $j=0; $j <= $#array ; $j++ )
     {
	    push (@string, split("\t",$array[$j]));
	    if ($column_number > $#string) {push (@column,undef);}
	    else {push (@column, $string[$column_number]);}
	    undef @string;
     }
    print STDERR join("\n", $column[0], $column[1], $column[2]), "\n";
    return (@column);
}

#-------------------------------------------------------------
sub calcMean {
    return sum(@_)/@_;
}

#-------------------------------------------------------------
sub calcStdDev {
	my(@array) = @_;
    my $n = $#array + 1;
	my $sum = sum(@array);
	$_ *= $_ for @array;
    my $sumOfSquares = sum(@array);
    my $variance = calcVariance( $sumOfSquares, $sum, $n );
    my $stddev = sqrt( $variance );
    return ($stddev,$variance);
}

#-------------------------------------------------------------
sub calcVariance {
	my ($sumOfSquares, $sum, $n) = @_;
	return ( $sumOfSquares - ( ( $sum * $sum ) / $n ) ) / ( $n - 1 );
}

#-------------------------------------------------------------

sub ReadFile {
    my ($in_file, $filename, $col_coords, $col_occup, $occupancy_hashref) = @_;
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
		or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
	}
	else { open( $inFH, "<", $in_file ) or die "error: $in_file cannot be opened:$!"; }

    my $buffer = "";
    my $sz_buffer = 0;
    my $timer2 = time();
    # counter for the markers we see
    my $marker_count = 0;
    
    my $regex_split_newline='\n';
    
    my $processed_memory_size = 0;
    my $offset=0;
    
    my @all_occups;
    
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
            my $pos = $string[$col_coords];
	    $pos+=0;
            my $occup = $string[$col_occup];
	    $occup+=0;
	    $occupancy_hashref->{$pos}->{$filename} =  $occup;
	    push(@all_occups,$occup);
	    undef @string;
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
    my $mean_occup=calcMean(@all_occups);
    return($mean_occup);
}

#--------------------------------------------------------------------------
sub clean {

my $text = shift;

$text =~ s/\r//g;
$text =~ s/\n//g;
return $text;
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
	if ( ! -d $wd ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input directory $wd: $!\n"
		);
	}
	#if ( -e $outfile ) {
	#	pod2usage(
	#		-exitval => 2,
	#		-verbose => 1,
	#		-message => "'$outfile' exists in target folder\n"
	#	);
	#}

}
