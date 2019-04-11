#!/usr/bin/perl

=head1 NAME

average_replicates.pl - Calculates the average occupancy profile based on several replicate occupancy profiles 

=head1 SYNOPSIS

perl -w average_replicates.pl --dir=<path to working dir> --output=<path to results file> --coordsCol=0 --occupCol=1 --pattern="occ.gz" --printData --sum [--help] 

 Required arguments:
    --dir | -i         path to directory with aggregate profiles
    --output | -out    output table file name
	
 Options:
   define column numbers in the input occupancy files (Nr. of the very first column is 0):
    --coordsCol | -cC  chromosome coordinates column Nr. (default: -cC 1)
    --occupCol | -oC   cumulative occupancy column Nr. (default: -oC 0)

   additional parameters
    --pattern | -p     occupancy profile file name extension template (default: occ.gz)
    --printData | -d   print all input occupancy columns to the output file
    --sum | -s         print column with sum of all occupancies for each nucleotide
	--list | -l        text file containing coma-separated list of all replicates (full path)
    --gzip | -z        compress the output
    --help | -h        Help
    
 Example usage:
    perl -w average_replicates.pl --dir=/mnt/hdd01/myWD --output=/mnt/hdd01/myWD/occup_tab.txt --pattern="occ.gz" --coordsCol=0 --occupCol=1 --printData
	
	OR
	
    perl -w average_replicates.pl -i /mnt/hdd01/myWD -out /mnt/hdd01/myWD/occup_tab.txt -p "occ.gz" -cC 0 -oC 1 -d
    
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

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

my (%occupancy,%NormFactors);

my $wd;
my $output;
my $filename_pattern='occ.gz';
my $list_file;

my $coordsCol=0;
my $occupCol=1;
my $addData;
my $printSum;

my $needsHelp;
my $useGZ;

my $options_okay = &Getopt::Long::GetOptions(
	'dir|i=s' => \$wd,
	'output|out=s'   => \$output,
	'pattern|p=s' => \$filename_pattern,
	'list|l=s' => \$list_file,

	'coordsCol|cC=s' => \$coordsCol,
	'occupCol|oC=s' => \$occupCol,
	
	'printData|d' => \$addData,
	'sum|s' => \$printSum,
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
	print STDERR "GZIP support enabled\n";
}
else {
	print STDERR "ZGIP support disabled\n";
}


if($list_file) { print STDERR "loading replicates list file: $list_file\n"; }

my $tm = localtime;
print STDERR "-----------------------\n",
join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",
join(":",$tm -> [2],$tm -> [1],$tm -> [0]),
"\n-----------------------\n";

#load list of replicated experiments
my (@names,@files,@dirs, @all_files);

if ($list_file) {
	open(LIST_FILE, "<$list_file") or die "can't read from file $list_file: $!";
	my @lines=<LIST_FILE>;
	close (LIST_FILE);
	my $pattern="[\\t\\s,;]";
	@all_files=split($pattern,join("",@lines));
	@all_files = grep /\S/, @all_files;
} else {
	opendir(DIR, "$wd") or die $!;
	@all_files = readdir(DIR);
	closedir(DIR);
}

foreach my $file (sort @all_files){
  if ($file =~ m/.*\.$filename_pattern$/){
	push(@files, $file);
	my $file_name = basename($file,  "\.$filename_pattern");
	push(@names, $file_name);
	}
}

for (my $i=0; $i<=$#files; $i++) {
	my $file_name = $names[$i];
	my $file = $files[$i];
	$NormFactors{$file_name} = ReadFile("$wd/$file", $file_name, $coordsCol, $occupCol, \%occupancy, @names);
}
print STDERR "calculating StDev, Variance, Sum and average.\nResults will be saved to $output\n";

# open pipe to Gzip or open text file for writing
  my ($gz_out_file,$out_file,$OUT_FHs);
  $out_file = $output;
	if ($useGZ) {
		$out_file =~ s/(.*)\.gz$/$1/;
		$gz_out_file = $out_file.".gz";
		$OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
	}
	else {
		open $OUT_FHs, '>', $output or die "Can't open $output for writing; $!\n";
	}


if ($addData) {
	if ($printSum) { print $OUT_FHs join("\t","position","Mean","Sum","stdev","Rel.Error", @names);   }
	else { print $OUT_FHs join("\t","position","Mean","stdev","Rel.Error", @names); }
} else {
	 print $OUT_FHs join("\t","position","Mean","stdev","rel.err."),"\n";
}

my $header=0; my $old_position;
my $total_counts = keys %occupancy;
print STDERR "processing $total_counts entries...\n";
my $work_progress_step = int($total_counts/10);
my $current_progress = $work_progress_step;
my $j=0;
for my $position ( sort {$a<=>$b} keys %occupancy) {
	if($current_progress == $j) {print STDERR "$current_progress from $total_counts...\n";$current_progress+=$work_progress_step;}
	if($j==0) { $old_position= $position; }
	if($old_position < $position-1 ) {
		my @hkeys = ($old_position+1..$position);
		$occupancy{@hkeys}{@names}=0
	}
    $j++;
    my @temp; my $coord;
    for my $name ( sort keys %{ $occupancy{$position} }) {
			if ( ($header==0) && ($addData) ){ print $OUT_FHs "\t$name"; }
			my $occup;
			if ($NormFactors{$name}==0) {
				print STDERR "index:$j\tname:$name\tposition:$position\tnorm factor:$NormFactors{$name}\toccupancy:$occupancy{$position}{$name}\tratio:NA\n";
			}
			else { 	$occup=$occupancy{$position}{$name}/$NormFactors{$name};	}
		
			push(@temp, $occup);
    }
    if ($header == 0) { print $OUT_FHs "\n"; }
    $header=1;
    my $Mean=calcMean(@temp);
		my $sum=sum(@temp);
    my ($stdev,$variance)=calcStdDev(@temp);
    my $rel_err;
    if ($Mean!=0) {  $rel_err=$stdev/$Mean; }
    else {$rel_err=0;}


    if ( ($addData) && ($printSum) ) { print $OUT_FHs join("\t",$position,$Mean,$sum,$stdev,$rel_err,@temp),"\n"; }
    elsif ( ($addData) && (!$printSum) ) { print $OUT_FHs join("\t",$position,$Mean,$stdev,$rel_err,@temp),"\n"; }
    elsif ( (!$addData) && ($printSum) ) { print $OUT_FHs join("\t",$position,$Mean,$sum,$stdev,$rel_err),"\n"; }
    else { print $OUT_FHs join("\t",$position,$Mean,$stdev,$rel_err),"\n";   }
	
	#print STDERR join(" | ",$position,$Mean,$sum,$stdev,$rel_err,@temp),"\n";
	
    undef @temp;
	$old_position = $position;
}
print STDERR "done!\n";
$tm = localtime;
print STDERR "-----------------------\n",join("-",$tm -> [3],1+ $tm -> [4],1900 + $tm -> [5])," ",join(":",$tm -> [2],$tm -> [1],$tm -> [0]),"\n-----------------------\n";
close($OUT_FHs) or die $!;
exit();



###################################

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
	if( (!$n) or ($n<=1) ) { return (0,0); }
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
    my ($in_file, $filename, $col_coords, $col_occup, $occupancy_hashref, @names) = @_;
    my $filesize = -s $in_file; #determine file size in bytes
    my $size_counter_step=int($filesize/100);
    $filesize = int($filesize/1048576); # filesize in megabytes

    print STDERR "Reading $in_file file of $filesize MBs. Please wait...\n";

    #read file with by 4kb chanks
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
					if ($occup==0) { undef @string; next; }
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
            print STDERR "."; $processed_memory_size=0;
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
	if ( ( ! -d $wd ) && ( ! $list_file ) ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find either input directory $wd or list file $list_file: $!\n"
		);
	}

}
