#!/usr/bin/perl

=head1 NAME

aggregate_profile.pl - Calculates aggregate profile of sequencing read density around genomic regions

=head1 SYNOPSIS

perl -w aggregate_profile.pl --input=<in.occ.gz> --regions=<annotations.txt> [--expression=<gene_expression.rpkm>] --aligned=<output.aligned.tab.gz> --average_aligned=<output.aggregare.txt> [ --path2log=<AggregateProfile.log> --region_start_column=<column Nr.> --region_end_column=<column Nr.> --strand_column=<column Nr.> --chromosome_col=<column Nr.> --GeneId_column=<column Nr.> --Expression_columnID=<column Nr.> --Methylation_columnID=<column Nr.> --Methylation_columnID2=<column Nr.> --upstream_delta=<column Nr.> --downstream_delta==<column Nr.> --upper_threshold=<column Nr.> --lower_threshold=<column Nr.> --Methylation_threshold=<value|range_start-range_end> --overlap=<length> --library_size=<Nr.> --remove_minus_strand | --ignore_strand | --fixed_strand=[plus|minus] --invert_strand --input_occ --score --save_aligned --Cut_tail --chromosome=chrN --AgregateProfile --GeneLengthNorm --LibsizeNorm --PerBaseNorm --useCentre --use_default --verbose --help ]

 Required arguments:
 
  input files:
    --input | -in             Path to input *.OCC or *.OCC.GZ file; In combination with --ProcessDir flag path to folder with occupancy files for individual chromosomes to process. File name should contain a chromosome ID in the very beginning of every file name: chr1.Sample1.w1.occ.gz; 
    --regions | -reg          Regions annotation (genes annotation) table, containing information about alignment regions starts and stops as well as chromosome and strand information for gene lists (optional)

  output files:
    --aligned | -al           File name template for output table containing all aligned regions with their occupancy values; In combination with --ProcessDir flag path to output folder. File names will be generated automatically based on input file names and run settings: old_file_name.delta_1000_1000.w1.aligned.txt.gz
    --average_aligned | -av   File name template for output file containing aggregate profile data; In combination with --ProcessDir flag path to output folder. File names will be generated automatically based on input file names and run settings: old_file_name.delta_1000_1000.w1.average_aligned.txt
    
 Options:
 
  optional files:
    --expression | -exp       RPM expression values and flags
    --path2log | -log         path to program log file (default: AggregateProfile.log in working directory)
    
  define column numbers in the input regions annotation file (Nr. of the very first column is 0):
    --GeneId_column | -idC         GeneId column Nr. (default: 0)
    --Methylation_columnID | -m1C  "Nr. of reads where the nucleotide was methylated (5mC cytosines)" column Nr. (default: 2)  
    --Methylation_columnID2 | -m2C "Total Nr. of reads including methylated (CpG)" column Nr. (default: 3)
    --chromosome_col | -chrC       chromosome column Nr. (default: 6)
    --strand_column | -strC        strand column Nr. (default: 7)
    --region_start_column | -sC    region_start column Nr. (default:8)
    --region_end_column | -eC      region_end column Nr. (default:9)
 
  expression related flags and parameters (Nr. of the very first column is 0):
    --Expression_flag_columnID | -expFC    Expression flag column Nr. in expression file
    --Expression_value_columnID | -expVC   Expression value column Nr. in expression file
    --Expression_flag | -eFlag             Expression flag: remove genes marked "excluded" (default value) from calculations
    
  methylation-related column IDs and flags:
    --Apply_methylation_filter | -useMeth   activate the methylation filtering mode - uses for analysis only ranges with methylation within the range defined by --Methylation_threshold parameter
    --Methylation_threshold | -mT           Set the threshold for the methylation value. Keep only sites with methylation within the range. Default range: 50-100%.
    
  analysis-related flags and parameters:
    --upstream_delta | -upD         number of base pairs upstream from aligned starts considered in calculations (Default: 100)
    --downstream_delta | -downD     number of base pairs downstream from aligned starts considered in calculations (Default: 1500)
    --upper_threshold | -upT        set upper occupancy threshold to remove mapping/sequencing artifacts (Default: 10000)
    --lower_threshold | -loT        set lower occupancy threshold to avoid ORFs with low coverage (Default: 0)
    --overlap | -ov                 remove overlapping ranges from analysis: region starts should be located at least --overlap bases from each other (Default: 100)
    --library_size | -lS            sequencing library size (to use with -LibsizeNorm and --score)
    --chromosome | -chr             limit analysis to regions derived from specified chromosome only
    
    --useCenter | -uC               Use middle of the region for alignment instead of start of the region
    --remove_minus_strand  | -noMS  remove all genes marked as "minus" strand
    --ignore_strand | -noS          ignore strand information (mark all as "plus" strand)
    --fixed_strand | -fixS          ignore strand information and assign selected
    --invert_strand | -invS         invert start and stop strands
    --input_occ | -inOCC            use occupancy file as an input (*.occ or *.occ.gz)
    --score                         calculate RPM value (requires --library_size)
    --save_aligned | -sA            save aligned matrix
    --Cut_tail | -noTail            do not use reads downstream from regions end within the downstream_delta
    --window | -w                   running window parameter - should be equal to one, provided in bed2occupancy_average.pl parameter. (Default: 100; smallest possible: 1)
	
  normalization options
    --AgregateProfile | -aggr       calculates aggregate profile representing the average occupancy
    --GeneLengthNorm | -glN         normalize each profile to the region length (gene length)
    --LibsizeNorm | -lsN            perform sequencing library size normalization
    --PerBaseNorm | -pbN            compensate for the transcript length difference
	
  additional parameters
    --ProcessDir | -dir             generate aggregate profiles for each OCC file in the directory specified by --input flag
    --force                         overwrite output file with the same file name 
    --gzip | -z                     compress the output
    --verbose                       display additional run info
    --help| -h                      display detailed help
 
 
 Example usage:
 
Example 1:
 Generate aggregate profile based on the occupancy data for chromosome 1 only, extracted from Gzip-compressed occupancy file "name_template.occ.gz". The resulting aggregate profile will be saved to chr1.name_template.aggregate.txt and the underlying aligned matrix to chr1.name_template.aligned.txt.gz Whole data set will be normalized to the sequencing library size; aligned data will be compensated for a differences in transcript length; All reads downstream from transcription termination site but within the range of default "downstream_delta" (1500) will be ignored.
	
aggregate_profile.pl --input=name_template.occ.gz --regions=genes_annotations.txt --aligned=chr1.name_template.aligned.txt.gz --average_aligned=chr1.name_template.aggregate.txt --LibsizeNorm --PerBaseNorm --Cut_tail --library_size=20000000 --chromosome=chr1
	
 OR
	
aggregate_profile.pl -in chr1.name_template.occ.gz -reg genes_annotations.txt -al chr1.name_template.aligned.txt.gz -av chr1.name_template.aggregate.txt --LibsizeNorm --PerBaseNorm --Cut_tail -lS 20000000 -chr chr1

Example 2:
 Generate aggregate profile based on the occupancy data for chromosome 1 only, extracted from Gzip-compressed occupancy file "name_template.occ.gz". The resulting aggregate profile will be saved to chr1.name_template.aggregate.txt. Aligned matrix will not be saved.
 Align average occupancy plots at the center of a reference regions. No strand information available. Regions coordinates are taken from the "regions_annotations.txt" file. All calculation will be done only for bases with methylation in the range from 65 to 95% (Nr. of reads where the nucleotide was methylated divided to the total Nr. of reads)
	
aggregate_profile.pl --input=name_template.occ.gz --regions=regions_annotations.txt --dont_save_aligned --Apply_methylation_filter --Methylation_threshold=65-95 --useCenter --ignore_strand --average_aligned=Meth65_95.name_template.aggregate.txt --LibsizeNorm --library_size=20000000
    
 OR
    
aggregate_profile.pl -in name_template.occ.gz -reg regions_annotations.txt -discA -useMeth -mT=65-95 -uC -noS -av Meth65_95.name_template.aggregate.txt --LibsizeNorm -lS 20000000

Example 3:
 Generate aggregate profile for every file in the input folder:
 
aggregate_profile.pl -reg regions_annotations.txt -idC 0 -chrC 4 -strC 7 -sC 8 -eC 9 -al /mnt/output/Aggregate_profiles -av /mnt/output/Aggregate_profiles -in /mnt/input/Sample1/OCC -upD 1000 -downD 1000 --gzip --force -w 1 --PerBaseNorm --ProcessDir
	
command will scan input folder /mnt/input/Sample1/OCC for all files like chr1.sample1.w1.occ.gz, chr2.sample1.w1.occ.gz e.t.c. and save resulting aggregate profiles in the output folder /mnt/output/Aggregate_profiles using default file names like 

=head1 DESCRIPTION

=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 bed2occupancy_average.pl

 bed2occupancy_average.pl takes as input a bed file with coordinates of genomic features (promoters, enhancers, chromatin domains, TF binding sites, etc), and the files with continuous chromosome-wide occupancy (nucleosome occupancy, TF distribution, etc). Calculates normalized occupancy profiles for each of the features, as well as the aggregate profile representing the average occupancy centerd at the middle of the feature
 
=head1 AUTHORS

=over

=item 
 Yevhen Vainshtein <yevhen.vainshtein@igb.fraunhofer.de>
 
=item 
 Vladimir Teif
 
=back

=head2 Last modified

 16 October 2016
 
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


use strict "vars";
use Config;
use Time::localtime;
use Time::Local;
#use Roman;
use List::Util qw(sum);
use List::Util qw(first);
use Getopt::Long;
use Pod::Usage;
use IO::Dir;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

#  Time count Initialisation
my $timer1 = time();
my $tm = localtime;
my $start_sec = $tm -> [0];
my $start_min = $tm ->[1];
my $start_hour = $tm ->[2];
my $start_time = time();

# Default parametrs
my $verbose;
my $delta_1 = 100;
my $delta_2 = 1500;
my $upper_occup_threshold = 1000000;
my $lower_occup_threshold = 0;
my $overlap = 0;
my $remove_minus_strand;
my $ignore_strand;
my $invert_strand;
my $normalize;
my $GeneLengthNorm;
my $library_size_normalization;
my $apply_DivSum_normalization;
my $PerBaseNorm;
my $Cut_tail;
my $save_aligned;
my $apply_methylation_filter;
my $use_center;
my $window = 100;

my $GeneId_column = 0;
my $chromosome_nr_col=6;
my $strand_column=7;
my $region_start_column=8;
my $region_end_column=9;

my $GeneExpression_flag_column = 1;
my $GeneExpression_value_column = 2;

my $methylation_column=2;
my $methylation_column2=3;
my $methylation_range_right=1000000000000000;
my $methylation_range_left=50;
my $Methylation_threshold;

my $in_file; 
my ($out_path1, $out_path2); 
my $fixed_strand="plus";
my $input_occ;
my $calc_score;
my $Chromosome;

my $Gene_annotation_table;
my $library_size;
my $expression_file;
my $Expression_flag="excluded";

my $run_log_path = "AggregateProfile.log";
my $needsHelp;
my $useGZ;
my $force_rewrite;
my $ProcessDir;


my $options_okay = &Getopt::Long::GetOptions(
	# input files
	'input|in=s' => \$in_file,
	'regions|reg=s'   => \$Gene_annotation_table,
	'expression|exp=s'   => \$expression_file,
	#output files
	'aligned|al=s' => \$out_path1,
	'average_aligned|av=s' => \$out_path2,
	'path2log|log' => \$run_log_path,
	# column IDs in annotation file
	'region_start_column|sC=s' => \$region_start_column,
	'region_end_column|eC=s'   => \$region_end_column,
	'strand_column|strC=s' => \$strand_column,
	'chromosome_col|chrC=s'   => \$chromosome_nr_col,
	'GeneId_column|idC=s'   => \$GeneId_column,
	# expression flag column ID in expression file
	'Expression_columnID|expFC=s'   => \$GeneExpression_flag_column,
	'Expression_value_columnID|expVC=s' => \$GeneExpression_value_column,
	'Expression_flag|eFlag=s'   => \$Expression_flag,
	# methilation-related column IDs and flags
	'Apply_methylation_filter|useMeth' => \$apply_methylation_filter,
	'Methylation_columnID|m1C=s'   => \$methylation_column,
	'Methylation_columnID2|m2C=s' => \$methylation_column2,
	'Methylation_threshold|mT=s'   => \$Methylation_threshold,
	# analysis settings
	'upstream_delta|upD=s' => \$delta_1,
	'downstream_delta|downD=s' => \$delta_2,
	'upper_threshold|upT=s' => \$upper_occup_threshold,
	'lower_threshold|loT=s' => \$lower_occup_threshold,
	'overlap|ov=s' => \$overlap,
	'library_size|lS=s' => \$library_size,
	'chromosome|chr=s' => \$Chromosome,
	#flags
	'useCenter|uC' => \$use_center,
	'remove_minus_strand|noMS' => \$remove_minus_strand, # remove all genes marked as "minus" strand
	'ignore_strand|noS' => \$ignore_strand, # ignore strand information (mark all as "plus" strand)
	'fixed_strand|fixS=s' => \$fixed_strand,  # ignore strand information and assign selected
	'invert_strand|invS' => \$invert_strand, # invert start and stop strands
	'input_occ|inOCC' => \$input_occ, # use occupancy file as an input (*.occ
	'score' => \$calc_score, # calculate RPM value 
	'save_aligned|sA' => \$save_aligned,
	'Cut_tail|noTail' => \$Cut_tail,
	#normalization options
	'AgregateProfile|aggr' => \$normalize,
	'GeneLengthNorm|glN' => \$GeneLengthNorm,
	'LibsizeNorm|lsN' => \$library_size_normalization,
	'PerBaseNorm|pbN' => \$PerBaseNorm,
	'verbose' => \$verbose,
	'gzip|z' => \$useGZ,
	'force' => \$force_rewrite,
	'window|w=s' => \$window,
	'ProcessDir|dir' => \$ProcessDir,

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
	print STDERR "ZGIP support enabled\n";
}
else {
	print STDERR "ZGIP support disabled\n";
	if ( ($in_file =~ (/.*\.gz$/))  and (!$useGZ) ) {
		print STDERR "======================================\n";
		print STDERR "WARNING! Input files are probably compressed!\n";
		print STDERR "Use --gzip parameter to enable support for file compression/decompression!";
		print STDERR "======================================\n";
		exit;
	}
}

# set flags
$input_occ = $input_occ ? "yes" : "no";
$Cut_tail = $Cut_tail ? "yes" : "no";

if (($apply_methylation_filter) && ($Methylation_threshold =~ /(.*)-(.*)/ ) ) { $methylation_range_left=$1; $methylation_range_right=$2; }
elsif (($apply_methylation_filter) && ($Methylation_threshold =~ /(.*)/ ) ) { $methylation_range_left=$1; }
	   
#read arguments from command line

if ((!$library_size) && ($library_size_normalization)) {
    #code
    warn "please specify library size!";
    exit;
}

if ((!$library_size) && ($calc_score)) {
    #code
    warn "please specify library size!";
    exit;
}

#$out_path1 = $out_path1.".delta_".$delta_1."_".$delta_2.".txt.gz";
#$out_path2 = $out_path2.".delta_".$delta_1."_".$delta_2.".txt";
# Display input parametrs
print STDERR "======================================\n";
print STDERR "----------- path to files ------------\n";
print STDERR "run log file: $run_log_path\n";
print STDERR "in file:",$in_file, "\n";
if ($ProcessDir) { print STDERR "generate aggregate profiles for each OCC file in the $in_file directory\n"; }
print STDERR "aligned file:",$out_path1, "\n";
print STDERR "average file:", $out_path2, "\n";
print STDERR "gene annotation file: ",$Gene_annotation_table, "\n";
if ( defined $expression_file) { print STDERR "expression values and flags file: ",$expression_file, "\n"; }
print STDERR "---- annotation file columns IDs -----\n";
print STDERR "region ID column: ",$GeneId_column, "\n";
print STDERR "expression column (log2ratio or absolute): ",$GeneExpression_flag_column, "\n";
print STDERR "region start column: ",$region_start_column, "\n";
print STDERR "region end column: ",$region_end_column, "\n";
print STDERR "strand column: ",$strand_column, "\n";
print STDERR "--------- analysis settings ----------\n";
if ($invert_strand) { print STDERR "invert strand\n"; }
print STDERR "chromosome column: ",$chromosome_nr_col, "\n";
if ($apply_methylation_filter) {
	print STDERR "apply methylation filter:\n";
	print STDERR "methylation columns: ",$methylation_column, ", ",$methylation_column2,"\n";
	print STDERR "methylation range: ",$methylation_range_left, " to ", $methylation_range_right,"\n";	
}
if ($Chromosome) { print STDERR "process chromosome $Chromosome only\n"; }
print STDERR "delta_minus: ",$delta_1, "\n";
print STDERR "delta_plus: ",$delta_2, "\n";
if ($window < 1) { print STDERR "running window should be at least 1 (single nucleotide resolution)\nExiting...\n"; exit;}
if ($window =~ /\D/) { print STDERR "running window should be an integer value bigger than 0\nExiting...\n"; exit;}
print STDERR "running window size: $window \n";
print STDERR "upper occupancy threshold: ",$upper_occup_threshold, "\n";
print STDERR "lower occupancy threshold: ",$lower_occup_threshold, "\n";
print STDERR "allowed regions overlap (bp): ",$overlap, "\n";
print STDERR "replace reads by 0 downstream from region end: $Cut_tail\n";
if ($use_center) { print STDERR "align regions at the center\n"; }
if ($remove_minus_strand) { print STDERR "remove transcripts on minus-strands\n"; }
if ($ignore_strand) { print STDERR "ignore strand information (assume all on $fixed_strand)\n"; }
if ($save_aligned) { print STDERR "Save aligned occupancy profiles to tab-delimited text file (GZIP compressed when possible)\n"; }
print STDERR "------- Normalization options: -------\n";
if ($normalize) { print STDERR "normalize occupancy to a number of regions starts\n"; }
if ($GeneLengthNorm) { print STDERR "normalize each region occupancy to the region length\n"; }
if ($library_size_normalization) { print STDERR "Apply library size normalization\n";}
if ($calc_score) { print STDERR "calculating RPM values\n"; }
if ($PerBaseNorm) { print STDERR "Normalize aggregate profile by regions Nr. at each base\n"; }
if ($verbose) { print STDERR "\n----------------------------\n print service information to the console\n"; }
print STDERR "======================================\n";

# exit script if output files exist
if (! $force_rewrite) {
    if ( -e $out_path1 ) {
	    print STDERR "output file $out_path1 is exists already! \nExiting\n";
	    exit;
    }
    elsif ( -e $out_path2 ) {
	    print STDERR "output file $out_path2 is exists already! \nExiting\n";
	    exit;
    }
}


my ($LIST_array_ref, $TS_positions_ref, $TE_positions_ref, $chromosomes_ref, $GeneIDs_ref, $annotatiom_ref) = readAnnotation($Gene_annotation_table,$chromosome_nr_col,$GeneId_column,$methylation_column,$region_start_column,$region_end_column,$methylation_column2,$strand_column,$Chromosome,$ignore_strand,$invert_strand);
# read columns with transcription start position, starnd, chromosomes

my @LIST_array  = @{$LIST_array_ref};
my @TS_positions = @{$TS_positions_ref};
my @TE_positions = @{$TE_positions_ref};
#my @chromosomes = @{$chromosomes_ref};
#my @GeneIDs = @{$GeneIDs_ref};
#my %annotation= %{$annotatiom_ref};


my (@Expression,@Methylation_col1,@Methylation_col2,@Methylation);

if ($apply_methylation_filter) {
    @Methylation_col1 = Read_column($methylation_column,"5mC",1,\@LIST_array); # methylated cytosines
    @Methylation_col2 = Read_column($methylation_column2,"CpG",1,\@LIST_array); # all cytosines
    @Methylation = map { $Methylation_col1[$_] / ( $Methylation_col2[$_] + 0.0001 ) } 0..$#Methylation_col1;
}


my @strands;


if(! $ignore_strand ) {
    @strands = Read_column($strand_column,"strand", 1, \@LIST_array);
    for (my $i=1; $i<=$#strands; $i++) {
	
	if ( ($strands[$i] eq "plus") or ($strands[$i] eq "+") or ($strands[$i] eq "1")){ $strands[$i] = "plus";	}
	elsif ( ($strands[$i] eq "minus") or  ($strands[$i] eq "-") or ($strands[$i] eq "-1")) { $strands[$i] = "minus";	}
	else { warn "unidentified strand for $i ( $LIST_array[$i] ): $strands[$i]\n";}
    }
}
else { @strands = ("$fixed_strand") x $#TS_positions; }

if ($invert_strand) {
    for (my $i=1; $i<=$#strands; $i++) {
	
	if ( ($strands[$i] eq "plus") or ($strands[$i] eq "+") ){ $strands[$i] = "minus";	}
	elsif ( ($strands[$i] eq "minus") or  ($strands[$i] eq "-") ) { $strands[$i] = "plus";	}
	else { warn "unidentified strand for $i ( $LIST_array[$i] ): $strands[$i]\n";}
    }
}

#------------------------------------------------------------------------------
# read expression flag and values from translatome output file (if specified)
my %expression_flags;
if ( defined $expression_file) {
    print STDERR "Reading expression flags from $expression_file file...\n";

    open(EXPRESSION_FLAGS, "$expression_file") or die "can't read from file $expression_file: $!";
    while (<EXPRESSION_FLAGS>) {
	for my $chank (split/\r\n/) {
	    my @words = split ("\t", clean($chank));
	    $expression_flags {$words[0]} = $words[$GeneExpression_flag_column];
	    $expression_flags {$words[0]}{'value'} = $words[$GeneExpression_value_column];
	    undef @words;
	    }		
    }
    close (EXPRESSION_FLAGS);
}


############################################

if (! $ProcessDir) {
	AGGREGATE( \%expression_flags, $in_file, $run_log_path, $out_path1, $out_path2, $delta_1, $delta_2, $verbose, $apply_methylation_filter, $remove_minus_strand, $use_center, $invert_strand, $Cut_tail, $normalize, $PerBaseNorm, $save_aligned, $Chromosome, $Gene_annotation_table,$chromosome_nr_col,$GeneId_column,$methylation_column,$region_start_column,$region_end_column,$methylation_column2,$strand_column ); }

elsif ($ProcessDir) {
    # process each *.bed file in the folder
    my (%dir, @dir_list, @text_list);
    my $start_dir = $in_file;
    tie %dir, IO::Dir, $start_dir;
    foreach (keys %dir) { push (@dir_list, $_); }
	
	print STDERR "following files from the $start_dir folder will be processed:\n";
	foreach my $file (@dir_list) {
		if ($file =~ m/.*\.occ$/) { print STDERR $file, "\n";			}
        elsif ($file =~ m/.*\.occ.gz$/) { print STDERR $file, "\n";		}
		else { 	next; }
	}
	
	my $out_dir1 = $out_path1;
	my $out_dir2 = $out_path2;
    
	foreach my $file (@dir_list) {
		next unless (-f "$start_dir/$file");
		if ($file =~ m/.*\.occ$/) {
            $out_path1 = $file; $out_path1 =~ s/(.*)\.occ$/$1\.delta\_$delta_1\_$delta_2\.w$window\.aligned\.txt/;
            $out_path2 = $file; $out_path2 =~ s/(.*)\.occ$/$1\.delta\_$delta_1\_$delta_2\.w$window\.average_aligned\.txt/;
			}
        elsif ($file =~ m/.*\.occ.gz$/) {
            $out_path1 = $file; $out_path1 =~ s/(.*)\.occ.gz$/$1\.delta\_$delta_1\_$delta_2\.w$window\.aligned\.txt.gz/;
            $out_path2 = $file; $out_path2 =~ s/(.*)\.occ.gz$/$1\.delta\_$delta_1\_$delta_2\.w$window\.average_aligned\.txt/;
		}
		else {
			next;
		}
		
		# extract chromosome from the file name
		$file =~ /.*chr(?:_){0,1}([\dMTXY]*).*/i;
		$Chromosome = $1;

		$out_path1 = $out_dir1."/".$out_path1;
		$out_path2 = $out_dir2."/".$out_path2;
		AGGREGATE(\%expression_flags, "$start_dir/$file", $run_log_path, $out_path1, $out_path2, $delta_1, $delta_2, $verbose, $apply_methylation_filter, $remove_minus_strand, $use_center, $invert_strand, $Cut_tail, $normalize, $PerBaseNorm, $save_aligned, $Chromosome, $Gene_annotation_table,$chromosome_nr_col,$GeneId_column,$methylation_column,$region_start_column,$region_end_column,$methylation_column2,$strand_column );         
    }
}




$tm = localtime;
my $stop_sec = $tm -> [0];
my $stop_min = $tm ->[1];
my $stop_hour = $tm ->[2];;
my $stop_time = time();
my $message;
my $duration = $stop_time-$start_time;

$message = "\nStarted:\t$start_hour:$start_min:$start_sec\nnow:\t$stop_hour:$stop_min:$stop_sec\nduration:\t$duration sec.\n";

print STDERR "$message\nJob finished!\nBye!\n\n";

exit;



#------------------------------------------------------------------------------
# read annotation file
sub readAnnotation {
	#my (@LIST_array, @array, %annotation);
	my ($Gene_annotation_table,$chromosome_nr_col,$GeneId_column,$methylation_column,$region_start_column,$region_end_column,$methylation_column2,$strand_column,$Chromosome,$ignore_strand,$invert_strand) = @_;
	my @LIST_array;
	my %annotation;
	my $Annotated_lines=1;
	
	print STDERR "Reading $Gene_annotation_table file...\n";
	open(LIST_FILE, "<$Gene_annotation_table") or die "can't read from file $Gene_annotation_table: $!";
	
	while (<LIST_FILE>) {
		for my $chank  (split/\r\n/) {
			my $text = clean($chank);
			$Annotated_lines++;
			
			my @temp = split /[\t\r\f\n\,]/, $text;
			$annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'start'}=$temp[$region_start_column];
			$annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'stop'}=$temp[$region_end_column];
			if ($apply_methylation_filter) {
				$annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'methylation'}=$temp[$methylation_column] / ( $temp[$methylation_column2] + 0.0001 );
			}
			if ( (! $ignore_strand ) && (! $invert_strand) ){
				if ( ($temp[$strand_column] eq "plus") or ($temp[$strand_column] eq "+") or ($temp[$strand_column] eq "1")){
					$annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'strand'} = "plus";
					}
				elsif ( ($temp[$strand_column] eq "minus") or ($temp[$strand_column] eq "-") or ($temp[$strand_column] eq "-1")){
					$annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'strand'} = "minus";
					}
				else {
					warn "can't read strand info for $temp[$GeneId_column]\n";
				}
			}
			elsif ( (! $ignore_strand ) && ($invert_strand) ) {
				if ( ($temp[$strand_column] eq "plus") or ($temp[$strand_column] eq "+") or ($temp[$strand_column] eq "1")){
					$annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'strand'} = "minus";
					}
				elsif ( ($temp[$strand_column] eq "minus") or ($temp[$strand_column] eq "-") or ($temp[$strand_column] eq "-1")){
					$annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'strand'} = "plus";
					}
				else {
					warn "can't read strand info for ",$annotation{"$temp[$GeneId_column]"}{'strand'},"\n";
				}
			} else { $annotation{$temp[$chromosome_nr_col]}{"$temp[$GeneId_column]"}{'strand'} = "$fixed_strand"; }
	
	
			if(!$Chromosome) {
				push(@LIST_array, $text);
			} else {
				my $list_chr=$temp[$chromosome_nr_col];
				my $selected_chr=$Chromosome;
				$list_chr =~ s/chr//i;
				$selected_chr =~ s/chr//i;
				if ( ($list_chr == $selected_chr ) or ($list_chr eq $selected_chr) ) {
					push(@LIST_array, $text);
				}
			}
			
			undef @temp;
	
		}
	}
	close (LIST_FILE);
	
	if($Chromosome) {
		my $loaded_transcripts= keys %{ $annotation{$Chromosome } };
		print STDERR "chromosome $Chromosome: $loaded_transcripts from $Annotated_lines annotated features loaded\n";
	} else {
		print STDERR "annotation for $Annotated_lines features loaded\n";
	}
	
	my @TS_positions = Read_column($region_start_column,"Region start column",1,\@LIST_array);
	my @TE_positions = Read_column($region_end_column,"Region end column",1,\@LIST_array);
	my @chromosomes = Read_column($chromosome_nr_col,"Chromosome",1,\@LIST_array);
	my @GeneIDs = Read_column($GeneId_column,"Region ID",1,\@LIST_array);

	return(\@LIST_array,\@TS_positions, \@TE_positions, \@chromosomes, \@GeneIDs, \%annotation);
}
#------------------------------------------------------------------------------
# read file with occupancies
sub openOCC {
	my ($in_file,$min_coord,$max_coord,$Chromosome,$annotation_hash_ref) = @_;
	my %annotation= %{$annotation_hash_ref};

	my $filesize = -s $in_file; #determine file size in bytes
	my $size_counter_step=int($filesize/100);
	$filesize = int($filesize/1048576); # filesize in megabytes
	
	print STDERR "loading reads for the coordinates range [$min_coord..$max_coord] from $in_file file of $filesize MBs. \nPlease wait...\n";
	
	#read file with by 4kb chanks
	#@coord_occ_array=();
	my $BUFFER_SIZE = 1024*128;
	
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
	
	my $regex_split_tab='.*\t(.*)';
	my $regex_split_newline='\n';
	# occupancy file with chromosome ID as the very first column: chromosome | coordinate | occupancy
	my $regexp_pattern1 = '(^chr\S{1,2})\s(\d*)\s(\d*)$'; 
	# usual (per-chromosome) occupancy file: coordinate | occupancy
	my $regexp_pattern2 = '(^\d*)\s(\d*(?:\.)?.*)$';
	my $processed_memory_size = 0;
	my $offset=0;
	
	my %occupancy;
	my $false_counter=0;
	my $incomplete_lines_counter=0;
	my $total_counter=0;
	my $loaded_lanes=0;
	my $last_line;
	my $total_bytes_loaded=0;
	
	while ((my $n = read($inFH, $buffer, $BUFFER_SIZE)) !=0) {
		if ($n >= $BUFFER_SIZE) {
		$buffer .= $inFH;
		}
		my @lines = split(/$regex_split_newline/o, $buffer);
		my $NrOfLines = $#lines;
		if ( $last_line =~ /(.*)IO\:\:.*/ ) { $lines[0] = $1.$lines[0]; $incomplete_lines_counter++; }
		else { $last_line=""; }
		# process each line in zone file
		foreach my $line (@lines) {
			$total_counter++;
			if ($line =~ /$regexp_pattern1/) {
				if (!$3) { next; }
				my $chrom=$1;
				my $pos = $2; $pos+=0;
				if( ($pos > $max_coord) || ($pos < $min_coord) ) { next; }
				$loaded_lanes++;
				my $occup_val = $3; $occup_val+=0;
				$occupancy{$chrom}{$pos}=$occup_val;
			}
			elsif ($line =~ /$regexp_pattern2/) {
				my $pos = $1; $pos+=0;
				if( !$Chromosome || ($Chromosome && ($pos > $max_coord) || ($pos < $min_coord)) ) { next; }
				$loaded_lanes++;
				my $occup_val = $2; $occup_val+=0;
				$occupancy{$Chromosome}{$pos}=$occup_val;
			}
			else {
				$false_counter++;
			}
		}
		$last_line = $lines[$NrOfLines];
		$processed_memory_size += $n;
		$total_bytes_loaded += $n;
		$offset += $n;
		if(int($processed_memory_size/1048576)>= $filesize/10) {
			print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.                         \r";
			$processed_memory_size=0;
		}
		undef @lines;
		$buffer = "";
	}
	
	my $NrLanes = keys %{ $annotation{$Chromosome} };

	my $duration = time()-$timer2;
	print STDERR " done in ", time()-$timer2, " seconds.\n";
	my $total_mbytes_loaded = sprintf "%.2f", $processed_memory_size/1048576;
	print STDERR $total_mbytes_loaded, " Mbs processed \n\n";
	print STDERR $false_counter+$incomplete_lines_counter," strings from $total_counter failed to load. $incomplete_lines_counter strings recovered\n",
	"$NrLanes regions with annotation loaded for chromosome $Chromosome\n",
	$total_counter-$loaded_lanes-$false_counter," coordinates out of range\n";
	#exit;
	close($inFH) or die $!;
	
	return(%occupancy);
}




#-------------------------------------------------------------------------------
# aggregate profile

sub AGGREGATE {
	my ($expression_hash_ref, $in_file, $run_log_path, $out_path1, $out_path2, $delta_1, $delta_2, $verbose, $apply_methylation_filter, $remove_minus_strand, $use_center, $invert_strand, $Cut_tail, $normalize, $PerBaseNorm, $save_aligned, $Chromosome,$Gene_annotation_table,$chromosome_nr_col,$GeneId_column,$methylation_column,$region_start_column,$region_end_column,$methylation_column2,$strand_column ) = @_;
	
	my %expression_flags= %{$expression_hash_ref};
	
	my @average_occ_freq_distr;
	my @output_array;
	my @output_array2;
	my @non_zero_counter = (0) x ($delta_1+$delta_2+1);
	my $genes_counter=0;
	#my @splice_array;
	
	my $old_id="GeneID";
	my $old_transcript_start=0;
	my $ignore_overlap = "yes";
	my $old_strand="plus";
	
	#initialize counters
	my $removed_minus=0;
	my $removed_overlap=0;
	my $removed_threshold=0;
	my $removed_by_chromosome=0;
	my $removed_byID=0;
	my $removed_byExprFlag=0;
	my $removed_by_methylation_filter=0;
	my $removed_outofrange=0;
	
	my @length_array;
	
	# $Chromosome value is not defined here - should be filtered from the file name, when scanning directory for files
	my ($LIST_array_ref, $TS_positions_ref, $TE_positions_ref, $chromosomes_ref, $GeneIDs_ref, $annotatiom_ref) = readAnnotation($Gene_annotation_table,$chromosome_nr_col,$GeneId_column,$methylation_column,$region_start_column,$region_end_column,$methylation_column2,$strand_column,$Chromosome,$ignore_strand,$invert_strand);
	# read columns with transcription start position, starnd, chromosomes
	
	#my @LIST_array  = @{$LIST_array_ref};
	my @TS_positions = @{$TS_positions_ref};
	my @TE_positions = @{$TE_positions_ref};
	my @chromosomes = @{$chromosomes_ref};
	#my @GeneIDs = @{$GeneIDs_ref};
	my %annotation= %{$annotatiom_ref};
	# find max and min coord
	my ($max_coord, $min_coord);
	$max_coord=max(@TS_positions,@TE_positions)+2*$delta_2+1;
	$min_coord=min(@TS_positions,@TE_positions)-2*$delta_1-1;
	
	if($min_coord<0) {$min_coord=0;}

	print STDERR "Reading nucleosome occupancy for chromosome $Chromosome from $in_file...\n";
	
	my %occupancy = openOCC($in_file,$min_coord,$max_coord,$Chromosome, \%annotation);
	foreach my $chrom (sort { $a<=>$b || $a cmp $b } keys %occupancy) {
		my @occ_array_plus;
		my @occ_array_minus;
		my @splice_array;
		print STDERR "processing occupancy for chromosome $chrom...\n";
	
		for my $pos  ( keys %{ $occupancy{$chrom} } ) {
				$occ_array_plus[$pos] = $occupancy{$chrom}{$pos};
				$occ_array_minus[$pos] = $occupancy{$chrom}{$pos};
		}
		#my $chr_size= keys %{ $occupancy{$chrom} };
		undef %occupancy;
		
		#print STDERR "\nchromosome $chrom of $chr_size bases sorted\nCalculating aggregate profile...";
		
		my $work_progress_step = int($#TS_positions/10);
		my $current_progress = $work_progress_step;
		
		my $loaded_transcripts= keys %{ $annotation{$chrom} };
		my $j=0;
	
		foreach my $gene_id ( keys %{ $annotation{$chrom} } ) {
			my ($position_start, $position_end, $max_occup_counts, $transcript_length);
	
			if ($verbose) { print STDERR "\>\> $j from $loaded_transcripts\t"; }
			elsif($current_progress == $j) {print STDERR ".";$current_progress+=$work_progress_step;}
		
			# methylation filter
			if ($apply_methylation_filter) {
				if(($annotation{$chrom}{$gene_id}{'methylation'} <= $methylation_range_left) || ($annotation{$chrom}{$gene_id}{'methylation'} >= $methylation_range_right)) { $removed_by_methylation_filter++; $j++; next; }
			}	
			# check for repetitive transcript IDs.
			if($old_id eq $gene_id) { $removed_byID++; $j++; next; }
			else { $old_id=$gene_id; }    
		
			# remove genes by expression flags
			if ((keys %expression_flags) && ($expression_flags {$gene_id} eq $Expression_flag )) { $removed_byExprFlag++; $j++; next; }
			 
			 # remove minus-strand
			if ($remove_minus_strand) {	
				if ($annotation{$chrom}{$gene_id}{'strand'} eq "minus") { $removed_minus++; $j++; next; }
				else { $position_start=$annotation{$chrom}{$gene_id}{'start'};
					  $position_end=$annotation{$chrom}{$gene_id}{'stop'};}
			}
			
			if (!$annotation{$chrom}{$gene_id}{'strand'}) {
				#code
				$annotation{$chrom}{$gene_id}{'strand'}="plus";
			}
			
			if($old_strand eq $annotation{$chrom}{$gene_id}{'strand'}) { $ignore_overlap="no"; }
			else { $ignore_overlap="yes"; $old_strand=$annotation{$chrom}{$gene_id}{'strand'}; }
		
			
			my $check_TSS_position;
			# get transcription start position for plus and minus strand
			if ($annotation{$chrom}{$gene_id}{'strand'} eq "minus") {
				if ($annotation{$chrom}{$gene_id}{'start'} > $annotation{$chrom}{$gene_id}{'stop'} ) {
					$position_start=$annotation{$chrom}{$gene_id}{'stop'};
					$position_end=$annotation{$chrom}{$gene_id}{'start'};
				} else {
					$position_start=$annotation{$chrom}{$gene_id}{'start'};
					$position_end=$annotation{$chrom}{$gene_id}{'stop'};
				}
				$check_TSS_position=$old_transcript_start-$overlap;
			}
			elsif ($annotation{$chrom}{$gene_id}{'strand'} eq "plus") {
				$position_start=$annotation{$chrom}{$gene_id}{'start'};
				$position_end=$annotation{$chrom}{$gene_id}{'stop'};
				$check_TSS_position=$old_transcript_start+$overlap;
			}
			$transcript_length=abs($position_start-$position_end);
			
			if($position_end >= $#occ_array_minus ) { $j++; $removed_outofrange++; next; }
	
			if ($overlap==0) { $ignore_overlap="yes"; }
				
			# remove string if new TSS is to close to old one and transcript are on different strands
			if ($ignore_overlap eq "no") {
				if(( $check_TSS_position > $position_start ) && ($annotation{$chrom}{$gene_id}{'strand'} eq "plus")) { $removed_overlap++; $ignore_overlap="yes"; next; $j++; }
				elsif(( $check_TSS_position < $position_start ) && ($annotation{$chrom}{$gene_id}{'strand'} eq "minus")) { $removed_overlap++; $ignore_overlap="yes"; next; $j++; }
			}
			else { $old_transcript_start=$position_start; $ignore_overlap="yes";}
					
			#initialize array of 0 of 2*$delta+1 size 
			# read part of the array @occ_array; leave 0 if empty
			my ($start_of_region_occ, $end_of_region_occ);
			my ($start_of_region_splice, $end_of_region_splice);
			my ($start_shift,$end_shift);
			my ($old_central_point,$new_central_point);
			
			if ($use_center) {
				#code
				$new_central_point=int(($annotation{$chrom}{$gene_id}{'stop'}+$annotation{$chrom}{$gene_id}{'start'})/2);
			
				if($j==0) {$old_central_point=$new_central_point;}
				if ($overlap==0) { $ignore_overlap="yes"; }
				else  { $ignore_overlap="no"; }
			
				if ($ignore_overlap eq "no") {
				if(($old_central_point-$overlap < $new_central_point) && ($old_central_point+$overlap > $new_central_point ) )
				{ $removed_overlap++; $ignore_overlap="yes"; next; }
				}
				else { $old_central_point=$new_central_point; $ignore_overlap="yes";}
			
				$start_of_region_occ = $new_central_point-$delta_1;
				$end_of_region_occ = $new_central_point+$delta_2+1;
				$start_of_region_splice=0;
				$end_of_region_splice=$delta_1 + $delta_2 +1;
		
				push(@splice_array, @occ_array_plus[$start_of_region_occ..$end_of_region_occ]);
	
			} else {
				if ($Cut_tail eq "yes" ) {
					if($annotation{$chrom}{$gene_id}{'strand'} eq "plus") {
						my $shift=min($position_start+$delta_2,$position_end);
						
						$start_of_region_occ = $position_start-$delta_1;
						$end_of_region_occ = $shift+1;
						
						#shift splice array start if TSS-delta1<0
						if ($start_of_region_occ<0) {
						@splice_array = 0 x abs($start_of_region_occ);
						$end_of_region_splice=$end_of_region_occ;
						$start_of_region_occ=0;
						}
						else {
						$end_of_region_splice=$end_of_region_occ-$start_of_region_occ;
						}
						
						$start_of_region_splice=0;
						
						push(@splice_array, @occ_array_plus[$start_of_region_occ..$end_of_region_occ]);
						
						if ($end_of_region_splice < $delta_1+$delta_2+1) {
						#code
						push(@splice_array,0) for($end_of_region_splice+1..$delta_1+$delta_2+1);
						}
						
						if ($invert_strand) {
							my @temp = @splice_array;
							@splice_array = reverse(@temp);
							undef @temp;
						}
						
					}
					elsif ($annotation{$chrom}{$gene_id}{'strand'} eq "minus") {
						# position START and END swapped!!
						#my $shift = min($position_end+$delta_2,$position_start);
						#$start_of_region_occ = $position_end-$delta_2;
						#$end_of_region_occ=$shift+1;

						$start_of_region_occ = max($position_end-$delta_2,$position_start);
						$end_of_region_occ=$position_end+$delta_1+1;

						#shift splice array start if TSS-delta1<0
						if ($start_of_region_occ<0) {
						@splice_array = 0 x abs($start_of_region_occ);
						$end_of_region_splice=$end_of_region_occ;
						$start_of_region_occ=0;
						}
						else {
						$end_of_region_splice=$end_of_region_occ-$start_of_region_occ;
						}
								
						push(@splice_array, @occ_array_minus[$start_of_region_occ..$end_of_region_occ]);
								
						if (!$invert_strand) {
							my @temp = @splice_array;
							@splice_array = reverse(@temp);
							if ($end_of_region_splice < $delta_1+$delta_2+1) {
								#code
								push(@splice_array,0) for($end_of_region_splice+1..$delta_1+$delta_2+1);
							}
						}
						else {		    
							if ($end_of_region_splice < $delta_1+$delta_2+1) {
								#code
								my @temp = @splice_array;
								@splice_array = reverse(@temp);
								undef @temp;
								push(@splice_array,0) for($end_of_region_splice+1..$delta_1+$delta_2+1);
								@temp = @splice_array;
								@splice_array = reverse(@temp);
								undef @temp;
							}
						}
					}
				} 
				elsif ($Cut_tail eq "no" ) {
					#$start_of_region_splice=0;
					#$end_of_region_splice=$delta_1+$delta_2+1;
					if ($annotation{$chrom}{$gene_id}{'strand'} eq "plus") {
						
						$start_of_region_occ = $position_start-$delta_1;
						$end_of_region_occ = $position_start+$delta_2+1;
						
						push(@splice_array, @occ_array_plus[$start_of_region_occ..$end_of_region_occ]);
						
					}
					elsif ($annotation{$chrom}{$gene_id}{'strand'} eq "minus") {
						$start_of_region_occ = $position_end-$delta_2;
						$end_of_region_occ= $position_end+$delta_1+1;
						push(@splice_array, @occ_array_minus[$start_of_region_occ..$end_of_region_occ]);
						
						my @temp_array = @splice_array;
						@splice_array = reverse(@temp_array);
						undef @temp_array;
					}
					
					#print STDERR join ("\t",$gene_id,$chrom, $start_of_region_occ,$end_of_region_occ,$end_of_region_occ-$start_of_region_occ,@splice_array[0..5]),"\n";

				}
			}

			# check for artefacts
			$max_occup_counts = max(@splice_array);
			if (($max_occup_counts > $upper_occup_threshold) or ($max_occup_counts < $lower_occup_threshold)){ $removed_threshold++; $j++; next; }
				
			# normalize
			if ($GeneLengthNorm) {
				#code
				my $total_reads_per_transcript;
				$total_reads_per_transcript += $_ for @splice_array;
				my $norm_factor=$total_reads_per_transcript/$transcript_length;
				if ((!$total_reads_per_transcript) or ($total_reads_per_transcript == 0)) { @splice_array = (0) x ($delta_1+$delta_2+1); }
				else {
				my @temp_array = map { $_ / $norm_factor } @splice_array;
				undef @splice_array;
				@splice_array=@temp_array;
				undef @temp_array;
				}
			}
			
			
			if($library_size_normalization) {
				my $norm_factor=$library_size/1000000;
				if ((!$norm_factor) or ($norm_factor == 0)) { @splice_array = (0) x ($delta_1+$delta_2+1); }
				else {
				my @temp_array = map { $_ / $norm_factor } @splice_array;
				undef @splice_array;
				@splice_array=@temp_array;
				undef @temp_array;
				}
			}
		
			
			#if($apply_DivSum_normalization eq "yes") {
			#	my $norm_factor += $_ for @splice_array;
			#	if ($norm_factor == 0) { @splice_array = (0) x ($delta_1+$delta_2+1); }
			#	else {
			#		my @temp_array = map { $_ / $norm_factor } @splice_array;
			#		undef @splice_array;
			#		@splice_array=@temp_array;
			#		undef @temp_array;
			#		}
			#	}
			
			my $temp_string;
			if ($calc_score) {
				my $RPM = sum(@splice_array)/(1000000*$library_size);
				$temp_string = join("\t", $gene_id, $transcript_length, $RPM);
			}
			else {
				if ($window > 1) {
					my @temp_array;
					for (my $k=$window; $k<=$#splice_array; $k+=$window) {
						my $sum=sum(@splice_array[$k-$window..$k]);
						my $average = $sum/$window;
						push (@temp_array, $average);
					}
					undef @splice_array;
					@splice_array = @temp_array;
					undef @temp_array;
				}
				$temp_string = join("\t", $gene_id, $transcript_length, @splice_array);
			}
		
			push (@output_array, $temp_string);
			
			$temp_string = join("\t", $gene_id, $chromosomes[$j], $annotation{$chrom}{$gene_id}{'strand'}, $start_of_region_occ, $end_of_region_occ, $end_of_region_occ-$start_of_region_occ , $annotation{$chrom}{$gene_id}{'start'}, $annotation{$chrom}{$gene_id}{'stop'}, $transcript_length, $#splice_array);
			push (@output_array2, $temp_string);

			#increment each element of positioning array
			for (my $i=0; $i<= $#splice_array; $i++) {
				$average_occ_freq_distr[$i] += $splice_array[$i];
				$non_zero_counter[$i]++;
			}
			$genes_counter++;
			push(@length_array, $transcript_length);
			
			if ($verbose) {
				print STDERR $temp_string,"\n";
			}
			$j++; 
			undef @splice_array;
		}
		undef @occ_array_plus;
		undef @occ_array_minus;
	} 	# end for-loop

	print STDERR "\n\n==================\n",$#output_array+1," genomic regions evaluated\n";
	
	my @results;

	# generate average occupancy distribution
	if ( $normalize ) {
		for (my $k=0;$k<=$#average_occ_freq_distr; $k++ ) {
			my $normalized_occup_summ = $average_occ_freq_distr[$k] / $non_zero_counter[$k];
			push(@results,$normalized_occup_summ);
		}
	}
	else { @results = @average_occ_freq_distr; }
	
	undef @average_occ_freq_distr;
	@average_occ_freq_distr = @results;
	undef @results;
	
	#normalize to transcript nr at each base
	#@TrasncriptsNr_perBase
	if ( $PerBaseNorm ) {
		my @TrasncriptsNr_perBase;
		for (my $i=0; $i<=$#average_occ_freq_distr; $i++) {
			my $norm_factor=checkarray(\@length_array,$i);
			push(@TrasncriptsNr_perBase,$norm_factor);
		}
		
		for (my $k=0;$k<=$#average_occ_freq_distr; $k++ ) {
			if ($verbose) {
			print STDERR "$k:\t$average_occ_freq_distr[$k]\t$TrasncriptsNr_perBase[$k]\t";
			}
			if ($TrasncriptsNr_perBase[$k] == 0) {
			#code
			$TrasncriptsNr_perBase[$k]=1;
			}
			
			my $normalized_occup_summ = $average_occ_freq_distr[$k] / $TrasncriptsNr_perBase[$k];
			if ($verbose) {
			print STDERR "$normalized_occup_summ\n";
			}
			push(@results,$normalized_occup_summ);
		}
	}
	else { @results = @average_occ_freq_distr; }
	
	
	# writing outputs
	if ($save_aligned) {
		@output_array = grep /\S/, @output_array;
		
		# open pipe to Gzip or open text file for writing
		my $out_file = $out_path1;
		my $OUT_FHs;
		if ($useGZ) {
			$out_file =~ s/(.*)\.gz$/$1/;
			my $gz_out_file = $out_file.".gz";
			print STDERR "\nsaving aligned occupancy matrix (",$#output_array," rows) to $gz_out_file\n";
			$OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
		}
		else {
			print STDERR "\nsaving aligned occupancy matrix (",$#output_array," rows)to $out_path1\n";
			open $OUT_FHs, '>', $out_path1 or die "Can't open $out_path1 for writing; $!\n";
		}
	
		print $OUT_FHs join("\n", @output_array),"\n";
		close ($OUT_FHs);
		
		$out_file =~ s/\.txt/\.tab/;
		open (OCCUP_File, ">$out_file") or die "can't open file $out_file for writting: $!";
		print OCCUP_File join("\t", "Gene", "Chromosomes", "strand", "start_of_region_occ", "end_of_region_occ", "region length" , "start", "stop", "transcript_length", "array length"),"\n";
		print OCCUP_File join("\n", @output_array2),"\n";
		close (OCCUP_File);
	
	}
	
	print STDERR "\nsaving aggregate profile matrix to $out_path2\n";
	open (AveragedFile, ">$out_path2") or die "can't open file $out_path2 for writting: $!";
	my $first_shift =  first { $results[$_] >0 } 0..$#results;
	
	# generate coordinates
	my @coords;
	for (my $i=-$delta_1;$i<=$delta_2;$i+=$window) { push (@coords, $i); }
	
	for (my $k=0; $k<=($delta_1+$delta_1+1); $k++) {
		print AveragedFile $coords[$k],"\t",$results[$k],"\n";
	}
	
	close(AveragedFile);
	
	print STDERR "done\n\n$#output_array transcription starts has been identified.\n\n",
	"Removed transcripts:\n",
	"minus strand: $removed_minus\n",
	"TSS overlap: $removed_overlap\n",
	"threshold: $removed_threshold\n",
	"methylation: $removed_by_methylation_filter\n",
	"identical (by ID): $removed_byID\n",
	"removed by expression flag (excluded after CW): $removed_byExprFlag\n",
	"$removed_by_chromosome from $#TS_positions transcripts removed (other chromosomes)\n\n";
	
	open(LOG, ">>$run_log_path") or die "can't upend to a log file: $run_log_path for writing: $!";
	print LOG "done\n\n$#output_array transcription starts has been identified.\n\n",
	"Removed transcripts:\n",
	"minus strand: $removed_minus\n",
	"TSS overlap: $removed_overlap\n",
	"threshold: $removed_threshold\n",
	"methylation: $removed_by_methylation_filter\n",
	"identical (by ID): $removed_byID\n",
	"removed by expression flag (excluded after CW): $removed_byExprFlag\n",
	"$removed_by_chromosome from $#TS_positions transcripts removed (other chromosomes)\n\n";
	close (LOG);
	
	
	#undef @LIST_array;
	undef @TS_positions;
	undef @TE_positions;
	undef @strands;
	undef @chromosomes;
	#undef @GeneIDs;
	undef @Expression;
	#undef @splice_array;
	undef @non_zero_counter;
	undef @Methylation_col1;
	undef @Methylation_col2;
	undef @Methylation;

}



#------------------ non-zero counts at position N --------------------
sub checkarray {
    my $arref = shift;
    my $N = shift;
    my $counter=0;
    foreach my $entry (@{$arref}){
       if( $entry >= $N ){ $counter++; }
    }
    return ($counter);
}


#------------------ Read specified column --------------------
sub Read_column {
    my ($column_number, $column_name, $start_row, $array_ref) = @_;
    my @array = @{$array_ref};

    my (@column, @string);
    # read column of interest to the memory
    for(my $j=$start_row; $j <= $#array ; $j++ )
     {
	    push (@string, split("[\t\r\f\n\,]",$array[$j]));
	    if ($column_number > $#string) {push (@column,undef);}
	    else {push (@column, $string[$column_number]);}
	    undef @string;
     }
    print STDERR join("\n", $column_name, "------------", $column[0], $column[1], $column[2]), "\n\n";
    return (@column);
}

#--------------------------------------------------------------------------
sub clean {

my $text = shift;

$text =~ s/\r//g;
$text =~ s/\n//g;
return $text;
}
#--------------------------------------------------------------------------
sub max {
  my $max = $_[0];
  for ( @_[ 1..$#_ ] ) {
    $_ = 0 if !$_;
    $max = $_ if $_ > $max; }
  return($max);
}
#--------------------------------------------------------------------------
sub min {
  my $min = $_[0];
  for ( @_[ 1..$#_ ] ) {
    $_ = 0 if !$_;
    $min = $_ if $_ < $min; }
  return($min);
}
#--------------------------------------------------------------------------
sub average{
        my($data) = @_;
        if (not @$data) {
                #die("Empty array!\n");
				return undef;
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
#--------------------------------------------------------------------------
# catch exceptions
sub try(&) { eval {$_[0]->()} };
sub catch(&) { $_[0]->($@) if $@ }

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
			-message => "\nError specifying options.\n"
		);
	}
	if ( !-e $in_file and !-d $in_file ) {
			pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "\nCannot find either input file or input directory $in_file: $!'\n"
		);
	}
	if ( !-e $Gene_annotation_table ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "\nCannot find regions annotation file: '$Gene_annotation_table!'\n"
		);
	}

}
