#!/usr/bin/perl

=head1 NAME

bed2occupancy_average.pl - Calculates genome-wide occupancy based on the bed file with sequencing reads 

=head1 SYNOPSIS

perl -w bed2occupancy_average.pl --input=<in.bed.gz> --output=<out.occ.gz> [--chromosome_col=<column Nr.> --start_col=<column Nr.> --end_col=<column Nr.> --strand_col=<column Nr.> --window=<running window size> --consider_strand --ConvertAllInDir --help]

 Required arguments:
    -in       input BED file or BED.GZ file or directory containing bed or bed.gz files (if option -dir is used)
    -out      output occupancy file (OCC.GZ)
	
 Options:
 
 define column numbers in the input BED file (Nr. of the very first column is 0):
    -s        read start column Nr. (default: -s 1)
    -e        read end column Nr. (default: -e 2)
    -str      strand column Nr. (default: -str 5)
    -chr      chromosome column Nr. (default: -chr 0)
	-w        running window size (default: -w 100). Set to 0 to calculate frequencies for each base.

    -dir      set flag to convert all BED files in the directory to OCC
    -use      consider strand when calculating occupancy
    -help -h  Help
    
 Example usage:
    bed2occupancy_average.pl --input=chr1.name_template.bed --output=chr1.name_template.occ --window=1000 --consider_strand
    bed2occupancy_average.pl --input=DIR_WITH_BED --window=1000 --consider_strand --ConvertAllInDir
    bed2occupancy_average.pl --input=chr1.name_template.bed.gz --window=0 --consider_strand
	
	OR
    
	bed2occupancy_average.pl -in chr1.name_template.bed -out chr1.name_template.occ -w 1000 -use
    bed2occupancy_average.pl -in DIR_WITH_BED -w 1000 -use -dir
    bed2occupancy_average.pl -in chr1.name_template.bed.gz -w 0 -use

 Note:
    default column numbers refer to a standard bed file format:
	
	chromosome | read start | read end | read ID | score | strand
	=============================================================
	
=head1 DESCRIPTION

=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 bed2occupancy_average.pl

 bed2occupancy_average.pl converts all or specified chromosome.bed files to reads occupancy files.
 The running window occupancy file (*.OCC) is a text file containing normalized reads frequency distribution
 along each chromosome for each running window. *.OCC file have the following format:
  
 relative chromosome coordinate | average occupancy in the running window
 ========================================================================
 
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
use IO::Dir;
use List::Util 'sum';
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError) ;

# Variables set in response to command line arguments
# (with defaults)

my ($infile,$outfile);

#flags
my $flag;
my $ConvertAllInDir;
my $needsHelp;

# columns
my $start_col=1;
my $end_col=2;
my $strand_col=5;
my $chromosome_col=0;

# parameters
my $running_window=100;
my $region_start=0;
my $region_end;


my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$infile,
	'output|out=s'   => \$outfile,
	
	'consider_strand|use' => \$flag,
	'ConvertAllInDir|dir' => \$ConvertAllInDir,
	
	'start_col|s=s' => \$start_col,
	'end_col|e=s'   => \$end_col,
	'strand_col|str=s' => \$strand_col,
	'chromosome_col|chr=s'   => \$chromosome_col,

	'window|w=s'   => \$running_window,
	'region_start|rs=s' => \$region_start,
	'region_end|re=s'   => \$region_end,
	
	'help|h'      => \$needsHelp
);

# set flags
$flag = $flag ? 1 : 0;
$ConvertAllInDir = $ConvertAllInDir ? 1 : 0;

# Check to make sure options are specified correctly and files exist
&check_opts();


# Display input parameters
print STDERR "======================================\n";
print STDERR "input file or folder: ",$infile, "\n";
print STDERR "output file (if applicable): ", $outfile, "\n";
print STDERR "================ flags ==============\n";
print STDERR "convert all files in the directory: ",$ConvertAllInDir, "\n";
print STDERR "use strand information while converting to occupancy: ",$flag, "\n";
print STDERR "running window: ",$running_window, "\n";
print STDERR "region start: ",$region_start, "\n";
print STDERR "region end: ",$region_end, "\n";
print STDERR "=== input bed files coumn numbers =====\n";
print STDERR "chromosome column: ",$chromosome_col, "\n";
print STDERR "read start column: ",$start_col, "\n";
print STDERR "read end column: ",$end_col, "\n";
print STDERR "strand column: ",$strand_col, "\n";
print STDERR "======================================\n";


if ($ConvertAllInDir == 0) {
	if(!$outfile) {
        $outfile = $infile;
        $outfile =~ s/(.*)\.bed(.gz)?$/$1\.w$running_window\.occ/;
	}
	Bed2Occup($infile, $outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag); } 
elsif ($ConvertAllInDir == 1) {
    # process each *.bed file in the folder
    my (%dir, @dir_list, @text_list);
    my $start_dir = $infile;
    tie %dir, IO::Dir, $start_dir;
    foreach (keys %dir) { push (@dir_list, $_); }

    foreach my $file (@dir_list) {
        if ($file =~ m/.*\.bed$/) {
            $outfile = $file;
            $outfile =~ s/(.*)\.bed$/$1\.w$running_window\.occ/;
            Bed2Occup($start_dir."/".$file, $start_dir."/".$outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag);           
        }
        if ($file =~ m/.*\.bed.gz$/) {
            $outfile = $file;
            $outfile =~ s/(.*)\.bed.gz$/$1\.w$running_window\.occ/;
            Bed2Occup($start_dir."/".$file, $start_dir."/".$outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag);           
        }
    }
}

exit;

#--------------------------------------------------------------------------

sub Bed2Occup {
    my ($infile_name, $outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag) = @_;
    
	my $infile;
	if ( $infile_name =~ (/.*\.gz$/) ) {
		$infile = IO::Uncompress::Gunzip->new( $infile_name )
		or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
	}
	else { open( $infile, "<", $infile_name ) or die "error: $infile_name cannot be opened:$!"; }
    
    my $buffer = "";
    my $sz_buffer = 0;
    my $timer2 = time();
    # counter for the markers we see
    my $marker_count = 0;
    
    my $regex_split_newline='\n';
    
    my $filesize_in_bytes = -s $infile_name; #determine file size in bytes
    my $size_counter_step=int($filesize_in_bytes/10);
    my $filesize = int($filesize_in_bytes/1048576); # filesize in megabytes
    
    print STDERR "Reading $infile_name file of $filesize MBs. Please wait...\n";
    my $processed_memory_size = 0;
    my $offset=0;
    my $not_zero_counter=0;
    my $string_counter=0;
    my $BUFFER_SIZE = 1024;
    my @occup=();
    my $old_coordinate=1;
        
    while ((my $n = read($infile, $buffer, $BUFFER_SIZE)) !=0) {
        if (($n >= $BUFFER_SIZE) or (($n == $filesize_in_bytes))) {
            $buffer .= <$infile>;
        }
        my @lines = split(/$regex_split_newline/o, $buffer);
        # process each line in zone file
        foreach my $line (@lines) {
            chomp($line);
            my @newline1=split(/\t/, $line);
            my $chr_name=$newline1[$chromosome_col];
            my $start_nuc=$newline1[$start_col];
            my $end_nuc=$newline1[$end_col];
            
            my $strand;
            if ($flag == 1) {
                $strand=$newline1[$strand_col];
                if ($strand eq "-") { ($start_nuc, $end_nuc) = ($end_nuc, $start_nuc); }
            }
            
            my $nuc_length=$end_nuc-$start_nuc;
            my $occup_counter = $occup[$end_nuc];
            
            if(!$occup_counter) {
                $occup[$end_nuc]=0;
            }
             for (my $j=$start_nuc; $j<=$end_nuc; $j++) {
                if(!$occup[$j]) {$occup[$j]=0;}
                $occup[$j]++;
            }
        }
        $processed_memory_size += $n;
        $offset += $n;
        if(int($processed_memory_size/1048576)>= $filesize/10) {
            print STDERR "."; $processed_memory_size=0;
            }
        undef @lines;
        $buffer = "";
    }
    close($infile);
    
    print STDERR "\nsaving occupancy data.\n";

	# open pipe to Gzip or open text file for writing
	my $out_file = $outfile;
	$out_file =~ s/(.*)\.gz$/$1/;
	my $gz_out_file = $out_file.".gz";
	my $OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";

    $timer2 = time();
 
    #initialize running average for first 2*$running_window+1 nucleotides
    my ($sum, $normalized_occupancy);
    
    # modify running average by shifting
    if (!$region_end) {	$region_end=$#occup; }
	if ($running_window > 0 ) {
    print STDERR "calculating and printing normalized occupancy with a running window +/- ",$running_window,"\nPlease wait...";
		for (my $i=$region_start; $i<$region_end; $i+=$running_window) {
			$sum=0;
			$sum=sum(@occup[$i..$i+$running_window]);
	
			my $average = $sum/$running_window;
			
			if ($average !=0) {
			#	print $OUT_FHs join ("\t", $i+$running_window,0),"\n";
			#}
			#else {
				$normalized_occupancy = $occup[$i+int($running_window/2)]/$average;
				print $OUT_FHs join ("\t", $i+$running_window,$average),"\n";	    
			}
		}
		print STDERR "done\n";
	}
	else {
		print STDERR "calculating and printing occupancy\nPlease wait...";
		for (my $i=$region_start; $i<$region_end; $i++) {			
			if ($occup[$i] !=0) {
			#	print $OUT_FHs join ("\t", $i+$running_window,0),"\n";
			#}
			#else {
				print $OUT_FHs join ("\t", $i,$occup[$i]),"\n";	    
			}
		}
		print STDERR "done\n";
	}


    close($OUT_FHs);    
    print STDERR "Job finished!\n$outfile generated in ", time()-$timer2, " seconds.\nGood bye!\n";

}


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
	if ( !-e $infile ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input BED file $infile: $!\n"
		);
	}
	if ( ( !$outfile ) and  ($ConvertAllInDir == 0) ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Please specify output occupancy file name!\n"
		);
	}

}
