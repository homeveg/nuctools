#!/usr/local/bin/perl

###==================================================================================================
### Calculates the average occupancy profile based on several replicate occupancy profiles
### 
### average_replicates.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================

use strict;
use Time::localtime;
use Time::Local;
use File::Basename;
use List::Util qw(sum);


my $usage = "$0 -input=\"path to working dir\" -output=\"path to results file\" -coordsCol=0 -occupCol=2 -printData\n";

my $wd;
my $coordsCol=0;
my $occupCol=1;
my %occupancy;
my %NormFactors;
my $output;
my $addData="no";

if (@ARGV != 0) {
    foreach my $comand_line_flag (@ARGV) {
	if ($comand_line_flag =~ /-input=(.*)/i) { $wd = $1; }
	if ($comand_line_flag =~ /-output=(.*)/i) { $output = $1; }
	if ($comand_line_flag =~ /-coordsCol=(.*)/i) { $coordsCol = $1; }
	if ($comand_line_flag =~ /-occupCol=(.*)/i) { $occupCol = $1; }
	if ($comand_line_flag =~ /-printData/i) { $addData = "yes"; }
    }
}
else {
    print STDERR $usage,"\n";
    exit;
}

my $tm = localtime;
print STDERR "-----------------------\n",join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),"\n-----------------------\n";

#check if folder exists
if (! -d $wd) {
    print STDERR "can't open $wd\n";
    print STDERR $usage,"\n";
    exit; 
}
else {
    opendir(DIR, "$wd") or die $!;
    my @all_files = readdir(DIR);
    closedir(DIR);
    my (@names,@files);
    my $filename_pattern=".*\.bed";
    
    foreach my $file (sort @all_files){
      if ($file =~ m/$filename_pattern/){
        push(@files, $file);
        my $filename = basename($file,  ".bed");
        push(@names, $filename);
	print STDERR "process $filename...\n";
	#    my ($in_file, $filename, $col_coords, $col_occup, $occupancy_hashref) = @_;
        $NormFactors{$filename} = ReadFile("$wd/$file", $filename, $coordsCol, $occupCol, \%occupancy);
        }
    }
    
    
}

print STDERR "calcualting StDev, Variance and average.\nResults will be saved to $output\n";
open(OUT,">$output") or die $!;
my $size = keys %occupancy;
print OUT join("\t","position","Mean","stdev","Rel.Error");
#print OUT join("\t","position","Rel.Error");
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
	if ( ($header==0) && ($addData eq "yes") ){ print OUT "\t$name"; }
	my $occup;
	if ($NormFactors{$name}==0) {
	    print STDERR "index:$j\tname:$name\tposition:$position\tnorm factor:$NormFactors{$name}\toccupancy:$occupancy{$position}{$name}\tratio:NA\n";
	}
	else {
	    $occup=$occupancy{$position}{$name}/$NormFactors{$name};

	}

	push(@temp, $occup);
    }
    if ($header == 0) { print OUT "\n"; }
    $header=1;
    my $Mean=calcMean(@temp);
    my ($stdev,$variance)=calcStdDev(@temp);
    my $rel_err;
    if ($Mean!=0) {  $rel_err=$stdev/$Mean; }
    else {$rel_err=0;}

    if ($addData eq "yes") { print OUT join("\t",$position,$Mean,$stdev,$rel_err,@temp),"\n";   }
    else {print OUT join("\t",$position,$Mean,$stdev,$rel_err),"\n";   }
    #if ($addData eq "yes") { print OUT join("\t",$position,$rel_err,@temp),"\n";   }
    #else {print OUT join("\t",$position,$rel_err),"\n";   }
    undef @temp;
}
print STDERR "done!\n";
$tm = localtime;
print STDERR "-----------------------\n",join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),"\n-----------------------\n";
close(OUT) or die $!;
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
    my $n = 0;
    my $sum = 0;
    my $sumOfSquares = 0;

    foreach my $x ( @_ ) {
        $sum += $x;
        $n++;
        $sumOfSquares += $x * $x;
    }

    #   Calculate the variance.
    my $variance = ( $sumOfSquares - ( ( $sum * $sum ) / $n ) ) / ( $n - 1 );
    my $stddev = sqrt( $variance );
    return ($stddev,$variance);
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
    
    # open original file
    open(INPUT, $in_file) or die "error: $in_file cannot be opened\n";
    my $buffer = "";
    my $sz_buffer = 0;
    my $timer2 = time();
    # counter for the markers we see
    my $marker_count = 0;
    
    my $regex_split_newline='\n';
    
    my $processed_memory_size = 0;
    my $offset=0;
    
    my @all_occups;
    
    while ((my $n = read(INPUT, $buffer, $BUFFER_SIZE)) !=0) {
        if ($n >= $BUFFER_SIZE) {
        $buffer .= <INPUT>;
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
            print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
            }
        undef @lines;
        $buffer = "";
    }
    
    my $duration = time()-$timer2;
    
    print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\ndone.\n";
    close(INPUT) or die $!;
    my $mean_occup=calcMean(@all_occups);
    return($mean_occup);
}
