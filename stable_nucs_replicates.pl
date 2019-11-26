#!/usr/bin/perl

=head1 NAME

stable_nucs_replicates.pl - Finds stable nucleosomes using all replicates for the same cell state

=head1 SYNOPSIS

perl -w stable_nucs_replicates.pl --input=<path to input DIR> --output=<out.bed> --chromosome=chr1 [-coordsCol=0 -occupCol=2 -StableThreshold=0.5 --printData ] [--help] 

 Required arguments:
    --inputDir  | -in    path to directory with occupancy files
    --outputS | -o1      output stable regions file
    --outputF | -o2      output fuzzy regions file
    --chromosome | -chr  chromosome ID
    --windowSize | -w    running window size. Use same value as for average occupancy calculation (default: 1)
    --MeanAbove | -ma    set a threshold on mean occupancy: discard fuzzy regions with occupancy below a threshold (default: 0 - save all) 

 Options:
 
   define column numbers in the input occupancy files (Nr. of the very first column is 0):
    --coordsCol | -cC   chromosome coordinates column Nr. (default: -cC 1)
    --occupCol  | -oC   cumulative occupancy column Nr. (default: -oC 0)
    
 additional parameters
    --printData  | -d       print all input occupancy columns to the output file
    --StableThreshold | -t  set threshold on relative error (St.Dev/Mean) to define stable nucleosomes: relative error below a threshold (default: 0.5)
    --FuzzyThreshold | -ft   set threshold on relative error (St.Dev/Mean) to define fuzzy nucleosomes: relative error above a threshold (default: as stable)
    --fileExtention | -p    input files extension (default: occ)
    --fuzzy | -f            save fuzzy regions (below threshold). If --FuzzyThreshold nor set explicitly, use a --StableThreshold value instead
	
    --gzip | -z             compress the output
    --help | -h             Help
	
 Example usage:
 
    perl -w stable_nucs_replicates.pl --input=/mnt/hdd01/myWD --output=out.bed.gz --chromosome=chr1	
    
=head1 DESCRIPTION

=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 stable_nucs_replicates.pl

 stable_nucs_replicates.pl Finds stable nucleosomes using all replicates for the same cell state
 
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

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }


my $wd;
my $coordsCol=0;
my $occupCol=1;
my %occupancy;
my $chr;
my $StableThreshold=0.5;
my $FuzzyThreshold;
my $windowSize=1;
my $MeanAbove=0;

my %NormFactors;
my $output1;
my $output2;
my $addData;
my $filename_pattern="occ";

my $useGZ;
my $needsHelp;
my $fuzzy;

my $options_okay = &Getopt::Long::GetOptions(
	'inputDir|in=s' => \$wd,
	'outputS|o1=s'   => \$output1,
	'outputF|o2=s'   => \$output2,
	'chromosome|chr=s' => \$chr,
	'fileExtention|p=s' => \$filename_pattern,

	'coordsCol|cC=s' => \$coordsCol,
	'occupCol|oC=s' => \$occupCol,
	'windowSize|w=s' => \$windowSize,
	'MeanAbove|ma=s' => \$MeanAbove,

	'printData|d' => \$addData,
	'StableThreshold|t=s' => \$StableThreshold,
	'FuzzyThreshold|ft=s' => \$FuzzyThreshold,
	'gzip|z' => \$useGZ,
	'fuzzy|f' => \$fuzzy,

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
	print STDERR "GZIP support enabled\n";
}
else {
	print STDERR "GZIP support disabled\n";
	if ( ($filename_pattern =~ (/.*\.gz$/))  and (!$useGZ) ) {
		print STDERR "======================================\n";
		print STDERR "WARNING! Input file probably compressed!\n";
		print STDERR "Use --gzip parameter to enable support for file compression/decompression!";
		print STDERR "======================================\n";
		exit;
	}

}

print STDERR "======================================\n";
print STDERR "input directory: ",$wd, "\n";
print STDERR "input files extensions: ",$filename_pattern, "\n";
print STDERR "output stable regions file: ",$output1, "\n";
print STDERR "output fuzzy regions file: ",$output2, "\n";
print STDERR "======================================\n";
print STDERR "chromosome name: ",$chr, "\n";
print STDERR "Occupancy column ID: ",$occupCol, "\n";
print STDERR "Coordinates column ID: ",$coordsCol, "\n";
print STDERR "bin size: ",$windowSize, "\n";
print STDERR "stable threshold: rel.error. below ",$StableThreshold, "\n";
if($fuzzy) {
	print STDERR "save fuzzy data (below threshold): enabled\n";
	if(!$FuzzyThreshold) { $FuzzyThreshold=$StableThreshold;}
	print STDERR "fuzzy threshold: rel.error. above ",$FuzzyThreshold, "\n";
} else {
print STDERR "save fuzzy data (below threshold): disabled\n";
}
if($addData) {
print STDERR "save individual occupancy values: enabled\n";
} else {
print STDERR "save individual occupancy values: disabled\n";
}
# check if threshold on mean set
if ($MeanAbove == 0) {
	print STDERR "discard regions with occupancy below threshold: disabled\n";
} else {
	print STDERR "discard regions with occupancy below $MeanAbove: enabled\n";
}
print STDERR "======================================\n";

my $tm = localtime;
print STDERR "-----------------------\n",join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),"\n-----------------------\n";

opendir(DIR, "$wd") or die $!;
my @all_files = readdir(DIR);
closedir(DIR);
my (@names,@files);

foreach my $file (sort @all_files){
  if ($file =~ m/.*\.$filename_pattern$/){
	push(@files, $file);
	my $filename = basename($file,  "\.$filename_pattern");
	push(@names, $filename);
	}
}

print STDERR "process ",$#names+1," files. Please wait...\n";
for (my $i=0; $i<=$#files; $i++) {
	my $filename = $names[$i];
	my $file = $files[$i];
	$NormFactors{$filename} = ReadFile("$wd/$file", $filename, $coordsCol, $occupCol, \%occupancy, @names);
}


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
if ($fuzzy) {
	$out_file = $output2;
	if ($useGZ) {
		$out_file =~ s/(.*)\.gz$/$1/;
		$gz_out_file = $out_file.".gz";
		$OUT2_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
	}
	else {
		open $OUT2_FHs, '>', $output2 or die "Can't open $output2 for writing; $!\n";
	}	
}

print STDERR "calcualting StDev, Variance and average.\nResults will be saved to $output1";
if ($fuzzy) { print STDERR " and to $output2\n"; }
else {  print STDERR "\n"; }
my $size = keys %occupancy;
#print OUT join("\t","chr", "start", "stop", "Mean","stdev","Rel.Error");
my $header=0;
my $total_counts = keys %occupancy;
print STDERR "processing $total_counts entries...\n";
my $work_progress_step = int($total_counts/10);
my $current_progress = $work_progress_step;
my $j=0;
my $discarded_counter=0;
for my $position ( sort {$a<=>$b} keys %occupancy) {
    if($current_progress == $j) {
		print STDERR "$current_progress from $total_counts..."," "x40,"\r";
		$current_progress+=$work_progress_step;}
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
	if ( ($fuzzy) && ($rel_err > $StableThreshold) && ($rel_err!=0) && ($rel_err <= $FuzzyThreshold)) {
		if ( ($MeanAbove>0) && ($Mean<=$MeanAbove) ) { $discarded_counter++; next;} 
		if ($addData) { print $OUT2_FHs join("\t",$chr, $position-$windowSize, $position,$Mean,$stdev,$rel_err,@temp),"\n";   }
    	else {print $OUT2_FHs join("\t",$chr, $position-$windowSize, $position,$Mean,$stdev,$rel_err,),"\n";   }
	}
    elsif (  ($rel_err < $StableThreshold) && ($rel_err!=0) ) { 
    	if ($addData) { print $OUT1_FHs join("\t",$chr, $position-$windowSize, $position,$Mean,$stdev,$rel_err,@temp),"\n";   }
    	else {print $OUT1_FHs join("\t",$chr, $position-$windowSize, $position,$Mean,$stdev,$rel_err,),"\n";   }
    }
    undef @temp;
}
print STDERR "done!\n\n";
if ($MeanAbove>0) {
	print STDERR "-----------------------\n";
	print STDERR "$discarded_counter fuzzy regions discarded: mean occupancy below threshold $MeanAbove\n";
}
$tm = localtime;
print STDERR "-----------------------\n",
join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),
"\n-----------------------\n";
close($OUT1_FHs) or die $!;
if ($fuzzy) { close($OUT2_FHs) or die $!; }

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
    my ($in_file, $filename, $col_coords, $col_occup, $occupancy_hashref, @names) = @_;
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
	}
	else { open( $inFH, "<", $in_file ) or die "error: $in_file cannot be opened:$!"; }

    my $buffer = "";
    my $sz_buffer = 0;
    my $timer2 = time();
    # counter for the markers we see
    my $counter = 0;
    
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
			foreach my $file (@names) {
					if(! exists ( $occupancy_hashref->{$pos}->{$file} ) ) { $occupancy_hashref->{$pos}->{$file}=0; }
			}
			$occupancy_hashref->{$pos}->{$filename} =  $occup;
			push(@all_occups,$occup);
			undef @string;
			$counter++;
        }
        $processed_memory_size += $n;
        $offset += $n;
        if(int($processed_memory_size/1048576)>= $filesize/10) {
            $processed_memory_size=0;
			print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds"," "x40,"\r";
            }
        undef @lines;
        $buffer = "";
    }
    
    my $duration = time()-$timer2;
    
    print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\ndone.\n";
	  print STDERR "$counter positions added\n";

    close($inFH) or die $!;
    my $mean_occup=calcMean(@all_occups);
		undef @all_occups;
    return($mean_occup);
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
	if ( !$output1 ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify stable regions output file name\n"
		);
	}
	if ( !$chr ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify chromosome name\n"
		);
	}

}