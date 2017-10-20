#!/usr/bin/perl

=head1 NAME

bed2occupancy_average.pl - Calculates genome-wide occupancy based on the bed file with sequencing reads 

=head1 SYNOPSIS

perl -w bed2occupancy_average.pl --input=<in.bed.gz> --output=<out.occ.gz> [--outdir=<DIR_WITH_OCC> --chromosome_col=<column Nr.> --start_col=<column Nr.> --end_col=<column Nr.> --strand_col=<column Nr.> --window=<running window size> --consider_strand --ConvertAllInDir --help]

 Required arguments:
    --input | -in        input BED file or BED.GZ file or directory containing bed or bed.gz files (if option -dir is used)
    --output | -out      output occupancy file (OCC.GZ)
	
 Options:
 
 define column numbers in the input BED file (Nr. of the very first column is 0):
    --start_col | -s            read start column Nr. (default: -s 1)
    --end_col | -e              read end column Nr. (default: -e 2)
    --strand_col | -str         strand column Nr. (default: -str 5)
    --chromosome_col | -chr     chromosome column Nr. (default: -chr 0)
    --window | -w               running window size (default: -w 100). Set to 1 to calculate frequencies for each base.

    --ConvertAllInDir | -dir    set flag to convert all BED files in the directory to OCC
    --outdir | -odir            path to output folder (save to input dir if not specified)
    --consider_strand | -use    consider strand when calculating occupancy

    --gzip | -z                 compress the output
    --help | -h                 Help
    
 Example usage:
    bed2occupancy_average.pl --input=chr1.name_template.bed --output=chr1.name_template.occ --window=1000 --consider_strand
    bed2occupancy_average.pl --input=DIR_WITH_BED --outdir=DIR_WITH_OCC --window=1000 --consider_strand --ConvertAllInDir
    bed2occupancy_average.pl --input=chr1.name_template.bed.gz --window=0 --consider_strand
	
	OR
    
    bed2occupancy_average.pl -in chr1.name_template.bed -out chr1.name_template.occ -w 1000 -use
    bed2occupancy_average.pl -in DIR_WITH_BED -odir DIR_WITH_OCC -w 1000 -use -dir
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
use File::Basename;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

# Variables set in response to command line arguments
# (with defaults)

my ($infile,$outfile);
my $outdir;

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
my $running_window=1;
my $region_start=0;
my $region_end;
my $useGZ;


my $options_okay = &Getopt::Long::GetOptions(
	'input|in=s' => \$infile,
	'output|out=s'   => \$outfile,
	'outdir|odir=s'  => \$outdir,
	
	'consider_strand|use' => \$flag,
	'ConvertAllInDir|dir' => \$ConvertAllInDir,
	
	'start_col|s=s' => \$start_col,
	'end_col|e=s'   => \$end_col,
	'strand_col|str=s' => \$strand_col,
	'chromosome_col|chr=s'   => \$chromosome_col,

	'window|w=s'   => \$running_window,
	'region_start|rs=s' => \$region_start,
	'region_end|re=s'   => \$region_end,
	'gzip|z' => \$useGZ,
	
	'help|h'      => \$needsHelp
);

# set flags
$flag = $flag ? 1 : 0;

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
}

# Display input parameters
print STDERR "======================================\n";
print STDERR "input file or folder: ",$infile, "\n";
print STDERR "output file (if applicable): ", $outfile, "\n";
if ($outdir) { print STDERR "output folder: ", $outdir, "\n"; }
else { $outdir = dirname($infile); }
print STDERR "================ flags ==============\n";
if ($ConvertAllInDir) { print STDERR "convert all *.BED files in the $infile directory to *.OCC \n"; }
if ($flag == 1) { print STDERR "use strand information while converting to the occupancy \n"; }
print STDERR "running window: ",$running_window, "\n";
print STDERR "region start: ",$region_start, "\n";
print STDERR "region end: ",$region_end, "\n";
print STDERR "=== input bed files coumn numbers =====\n";
print STDERR "chromosome column: ",$chromosome_col, "\n";
print STDERR "read start column: ",$start_col, "\n";
print STDERR "read end column: ",$end_col, "\n";
print STDERR "strand column: ",$strand_col, "\n";
print STDERR "======================================\n";


if (! $ConvertAllInDir) {
	if(!$outfile) {
        $outfile = $infile;
        $outfile =~ s/(.*)\.bed(.gz)?$/$1\.w$running_window\.occ/;
	}
	BED_2_OCC($infile, $outdir."/".$outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag, $region_start, $region_end); } 
elsif ($ConvertAllInDir) {
    # process each *.bed file in the folder
    my (%dir, @dir_list, @text_list);
    my $start_dir = $infile;
    tie %dir, IO::Dir, $start_dir;
    foreach (keys %dir) { push (@dir_list, $_); }

    foreach my $file (@dir_list) {
        if ($file =~ m/.*\.bed$/) {
            $outfile = $file;
            $outfile =~ s/(.*)\.bed$/$1\.w$running_window\.occ/;
            BED_2_OCC($start_dir."/".$file, $outdir."/".$outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag, $region_start, $region_end);           
        }
        if ($file =~ m/.*\.bed.gz$/) {
            $outfile = $file;
            $outfile =~ s/(.*)\.bed.gz$/$1\.w$running_window\.occ/;
            BED_2_OCC($start_dir."/".$file, $outdir."/".$outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag, $region_start, $region_end);           
        }
    }
}

exit;

#--------------------------------------------------------------------------

sub BED_2_OCC {
    my ($infile_name, $outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag, $region_start, $region_end) = @_;
    
	my $infile;
	if ( $infile_name =~ (/.*\.gz$/) ) {
		$infile = IO::Uncompress::Gunzip->new( $infile_name )
		or die "IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip::GunzipError\n";
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
    my $line_counter=0;
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
			$line_counter++;
        }
        $processed_memory_size += $n;
        $offset += $n;
		
        if(int($processed_memory_size/1048576)>= $filesize/10) {
            print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.            \r";
            $processed_memory_size=0;
        }
        undef @lines;
        $buffer = "";
    }
    close($infile);
    
    print STDERR "\nsaving occupancy data.\n";

	# open pipe to Gzip or open text file for writing
	my $out_file = $outfile;
	my $OUT_FHs;
	
	if ($useGZ) {
		$out_file =~ s/(.*)\.gz$/$1/;
		my $gz_out_file = $out_file.".gz";
		$OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
	}
	else {
		open $OUT_FHs, '>', $outfile or die "Can't open $outfile for writing; $!\n";
	}

    $timer2 = time();
 
    # modify running average by shifting
    if (!$region_end) {	$region_end=$#occup; }
	my $LibSize_norm_factor = ($line_counter)/($region_end-$region_start);
    
    my $counter=0;
	if ($running_window > 1 ) {
		print STDERR "calculating and printing normalized occupancy with a running window +/- ",$running_window,"\nPlease wait...";
		for (my $i=$region_start; $i<$region_end; $i+=$running_window) {
			my $sum=sum(@occup[$i..$i+$running_window]);
			my $average = $sum/$running_window;
			my $normalized_occupancy = $average/$LibSize_norm_factor;
			print $OUT_FHs join ("\t", $i+$running_window,$normalized_occupancy),"\n";
            $counter++;
            
            if( $counter >= $line_counter/20 ) {
                print STDERR "."; $counter=0;
            }

		}
		print STDERR "done\n";
	}
	else {
		print STDERR "calculating and printing occupancy\nPlease wait...";
		for (my $i=$region_start; $i<$region_end; $i++) {			
			if ($occup[$i] !=0) {
				my $normalized_occupancy = $occup[$i]/$LibSize_norm_factor;
				print $OUT_FHs join ("\t", $i, $normalized_occupancy),"\n";
			}
            $counter++;
            if( $counter >= $line_counter/20 ) {
                print STDERR "."; $counter=0;
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
