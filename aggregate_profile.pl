#!/usr/local/bin/perl
###
###==================================================================================================
### Calculates aggregate profile of sequencing read density around genomic regions
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### aggregate_profile.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 19 July 2016
###==================================================================================================


use strict "vars";
use Config;
use Time::localtime;
use Time::Local;
#use Roman;
use List::Util qw(sum);

#  Time count Initialisation
my $timer1 = time();
my $tm = localtime;
my $start_sec = $tm -> [0];
my $start_min = $tm ->[1];
my $start_hour = $tm ->[2];
my $start_time = time();

# Default parametrs
my $verbose=0;
my $delta_1 = 100;
my $delta_2 = 1500;
my $upper_occup_threshold = 10000;
my $lower_occup_threshold = 0;
my $overlap = 0;
my $remove_minus_strand = "no";
my $ignore_strand="no";
my $invert_strand = "no";
my $normalize="no";
my $GeneLengthNorm ="no";
my $library_size_normalization = "no";
my $apply_DivSum_normalization="no";
my $PerBaseNorm = "no";
my $Cut_tail="no";
my $dont_save_aligned="no";
my $apply_methylation_filter = "no";
my $use_centre = "no";

my $chromosome_nr_col=6;
my $strand_column=7;
my $region_start_column=8;
my $region_end_column=9;
my $GeneId_column = 0;
my $GeneExpression_column ;

my $methylation_column=2;
my $methylation_column2=3;
my $methylation_range_right=1000000000000000;
my $methylation_range_left=50;

my $in_file; 
my ($out_path1, $out_path2); 
my $fixed_strand="plus";
my $input_occ="no";
my $calc_score="no";
my $Chromosome="chr1";

my $Gene_annotation_table;
my $library_size;
my $expression_file;

my $run_log_path = "AggregateProfile.log";

#read arguments from command line
if (@ARGV != 0) {
    foreach my $comand_line_flag (@ARGV) {
	#input files
	if ($comand_line_flag =~ /-input=(.*)/i) { $in_file = $1; } #raw reads, cental weigthed reads or occupancy data
        if ($comand_line_flag =~ /-regions=(.*)/i) { $Gene_annotation_table = $1; } # genes annotation table
        if ($comand_line_flag =~ /-expression=(.*)/i) { $expression_file = $1; } # RPM expreession values and flags (optional)
	
	#output files
	if ($comand_line_flag =~ /-aligned=(.*)/i) { $out_path1 = $1; }
	if ($comand_line_flag =~ /-average_aligned=(.*)/i) { $out_path2 = $1; }
        if ($comand_line_flag =~ /-path2log=(.*)/i) { $run_log_path = $1; }
	
	# column IDs in annotation file
        if ($comand_line_flag =~ /-region_start_column=(.*)/i) { $region_start_column = $1; }
        if ($comand_line_flag =~ /-region_end_column=(.*)/i) { $region_end_column = $1; }
        if ($comand_line_flag =~ /-strand_column=(.*)/i) { $strand_column = $1; }
        if ($comand_line_flag =~ /-chromosome_col=(.*)/i) { $chromosome_nr_col = $1; }
        if ($comand_line_flag =~ /-GeneId_column=(.*)/i) { $GeneId_column = $1; }
        if ($comand_line_flag =~ /-Expression_columnID=(.*)/i) { $GeneExpression_column = $1; }
	
	# methilation-related flags
	if ($comand_line_flag =~ /-Apply_methylation_filter/i) { $apply_methylation_filter="yes";}   
	if ($comand_line_flag =~ /-Methylation_columnID=(.*)/i) { $methylation_column = $1; }
	if ($comand_line_flag =~ /-Methylation_columnID2=(.*)/i) { $methylation_column2 = $1; }
	if ($comand_line_flag =~ /-Methylation_threshold=(.*)-(.*)/i) { $methylation_range_left=$1; $methylation_range_right=$2; }
	elsif ($comand_line_flag =~ /-Methylation_threshold=(.*)/i) { $methylation_range_left=$1; }
	
	# analysis settings
        if ($comand_line_flag =~ /-upstream_delta=(.*)/i) { $delta_1 = $1; } # sequence region before aligned starts
        if ($comand_line_flag =~ /-downstream_delta=(.*)/i) { $delta_2 = $1; } # sequence region after aligned starts
	if ($comand_line_flag =~ /-upper_threshold=(.*)/i) { $upper_occup_threshold = $1; } # set a threshold for maximum occupancy
	if ($comand_line_flag =~ /-lower_threshold=(.*)/i) { $lower_occup_threshold = $1; } # set a threshold for minimum occupancy
	if ($comand_line_flag =~ /-overlap=(.*)/i) { $overlap = $1; } # remove transcripts if starts closer than 100 (default) bases
	if ($comand_line_flag =~ /-library_size=(.*)/i) { $library_size=$1;} # library size 
	if ($comand_line_flag =~ /-chromosome=(.*)/i) { $Chromosome = $1; }
	if ($comand_line_flag =~ /-useCentre/i) { $use_centre = "yes"; }
		
	#flags
	if ($comand_line_flag =~ /-remove_minus_strand/i) { $remove_minus_strand="yes";} # remove all genes marked as "minus" starnd
	if ($comand_line_flag =~ /-ignore_strand/i) { $ignore_strand="yes";} # ignore starnd information (mark all as "plus" strand)
	if ($comand_line_flag =~ /-fixed_strand=(.*)/i) { $fixed_strand=$1;} # ignore starnd information and assign selected
	if ($comand_line_flag =~ /-invert_strand/i) { $invert_strand="yes";} # invert start and stop strands
	if ($comand_line_flag =~ /-input_occ/i) { $input_occ="yes";} # use occupancy file as an input (*.occ)
	if ($comand_line_flag =~ /-score/i) { $calc_score="yes";} # calculate RPM value 
	if ($comand_line_flag =~ /-dont_save_aligned/i) { $dont_save_aligned="yes";}
	if ($comand_line_flag =~ /-Cut_tail/i) { $Cut_tail="yes";}

	#normalization options
	if ($comand_line_flag =~ /-AgregateProfile/i) { $normalize="yes";}	
	if ($comand_line_flag =~ /-GeneLengthNorm/i) { $GeneLengthNorm="yes";}	
	if ($comand_line_flag =~ /-LibsizeNorm/i) { $library_size_normalization="yes";}
	if ($comand_line_flag =~ /-PerBaseNorm/i) { $PerBaseNorm="yes";}
	
        if ($comand_line_flag =~ /-use_default/i) { print STDERR "using default values: \n";}
        if ($comand_line_flag =~ /--verbose/i) { $verbose=1;}
	
	if ($comand_line_flag =~ /--help/i) {
	    print STDOUT <DATA>;
	    print STDOUT "\nPress <ENTER> button to exit... ";
	    <STDIN>;
	    exit;
	}
    }
}
else { warn
  "perl -w aggregate_profile.pl -input= -regions= -expression= ",
  "-aligned=[N] -average_aligned=[N] -path2log= ",
  "-region_start_column=[N] -region_end_column=[N] -strand_column=[N] -chromosome_col=[N] -GeneId_column=[N] -Expression_columnID=[N] ",
  "-Methylation_columnID=[N] -Methylation_columnID2=[N] ",
  "-upstream_delta=[N] -downstream_delta=[N] -upper_threshold=[N] -lower_threshold=[N] -Methylation_threshold=[N|N-N] -overlap=[N] -library_size=[N] ",
  "-remove_minus_strand | -ignore_strand | -fixed_strand=[plus|minus] -invert_strand -input_occ -score -dont_save_aligned -Cut_tail -chromosome=chrN ", 
  "-AgregateProfile -GeneLengthNorm -LibsizeNorm -PerBaseNorm -useCentre ",
  "-use_default --verbose --help    \n\n";
      exit;}

if ((!$library_size) && ($library_size_normalization eq "yes")) {
    #code
    warn "please specify library size!";
    exit;
}


$out_path1 = $out_path1.".delta_".$delta_1."_".$delta_2.".txt";
$out_path2 = $out_path2.".delta_".$delta_1."_".$delta_2.".txt";
# Display input parametrs
print STDERR "======================================\n";
print STDERR "path to a log file: $run_log_path\n";
print STDERR "======================================\n";
print STDERR "in file:",$in_file, "\n";
print STDERR "aligned file:",$out_path1, "\n";
print STDERR "average file:", $out_path2, "\n";
print STDERR "gene annotation file: ",$Gene_annotation_table, "\n";
if ( defined $expression_file) { print STDERR "expression values and flags file: ",$expression_file, "\n"; }
print STDERR "region ID column: ",$GeneId_column, "\n";
print STDERR "expression column (log2ratio or absolute): ",$GeneExpression_column, "\n";
print STDERR "region start column: ",$region_start_column, "\n";
print STDERR "region end column: ",$region_end_column, "\n";
print STDERR "strand column: ",$strand_column, "\n";
print STDERR "invert strand: ",$invert_strand, "\n";
print STDERR "chromosome column: ",$chromosome_nr_col, "\n";
print STDERR "methylation columns: ",$methylation_column, ", ",$methylation_column2,"\n";

print STDERR "apply methylation filter: ",$apply_methylation_filter, "\n";
print STDERR "methylation range: ",$methylation_range_left, " to ", $methylation_range_right,"\n";

print STDERR "selected chromosome: ",$Chromosome, "\n";
print STDERR "delta_minus: ",$delta_1, "\n";
print STDERR "delta_plus: ",$delta_2, "\n";
print STDERR "upper occupancy threshold: ",$upper_occup_threshold, "\n";
print STDERR "lower occupancy threshold: ",$lower_occup_threshold, "\n";
print STDERR "allowed regions overlap (bp): ",$overlap, "\n";
print STDERR "remove transcripts on minus-strands: $remove_minus_strand\n";
print STDERR "ignore strand information (assume all on $fixed_strand): $ignore_strand\n";
print STDERR "normalized occupancy to a number of TSS: $normalize\n";
print STDERR "Do not save aligned occupancy profiles: $dont_save_aligned\n";
print STDERR "Apply library size normalzation: $library_size_normalization\n";
print STDERR "Normalize by transcripts Nr. per base: $PerBaseNorm\n";
print STDERR "replace reads by 0 downstream from region end: $Cut_tail\n";
print STDERR "align regions at the centre: $use_centre\n";

print STDERR "======================================\n";

# read annotation file top
my (@LIST_array,
    #@coords_array,
    #@coord_occ_array,
    @array);
print STDERR "Reading $Gene_annotation_table file...\n";
open(LIST_FILE, "<$Gene_annotation_table") or die "can't read from file $Gene_annotation_table: $!";
while (<LIST_FILE>) { for my $chank  (split/\r\n/) { my $text = clean($chank); push(@LIST_array, $text); } }
close (LIST_FILE);

# read columns with transcription start position, starnd, chromosomes

my (@Expression,@Methylation_col1,@Methylation_col2,@Methylation);
my @TS_positions = Read_column($region_start_column,"Region start column",\@LIST_array);
my @TE_positions = Read_column($region_end_column,"Region end column",\@LIST_array);
my @chromosomes = Read_column($chromosome_nr_col,"Chromosome",\@LIST_array);
my @GeneIDs = Read_column($GeneId_column,"Region ID",\@LIST_array);

#@Expression = Read_column($GeneExpression_column,\@LIST_array);

if ($apply_methylation_filter eq "yes") {
    @Methylation_col1 = Read_column($methylation_column,"5mC",\@LIST_array); # methylated cytosines
    @Methylation_col2 = Read_column($methylation_column2,"CpG",\@LIST_array); # all cytosines
    @Methylation = map { $Methylation_col1[$_] / ( $Methylation_col2[$_] + 0.0001 ) } 0..$#Methylation_col1;
}


my @strands;


if($ignore_strand eq "no") {
    @strands = Read_column($strand_column,\@LIST_array);
}
else { @strands = ("$fixed_strand") x $#TS_positions; }

if ($invert_strand eq "yes") {
    for (my $i=1; $i<=$#strands; $i++) {
	
	if ($strands[$i] eq "plus") { $strands[$i] = "minus";	}
	elsif ($strands[$i] eq "minus") { $strands[$i] = "plus";	}
	else { warn "unidentifyed strand for $i ( $LIST_array[$i] ): $strands[$i]\n";}
    }
}

#------------------------------------------------------------------------------
# read expression flag from translatome output file (if specified)
my %expression_flags;
if ( defined $expression_file) {
    print STDERR "Reading expression flags from $expression_file file...\n";

    open(EXPRESSION_FLAGS, "$expression_file") or die "can't read from file $expression_file: $!";
    while (<EXPRESSION_FLAGS>) {
	for my $chank (split/\r\n/) {
	    my @words = split ("\t", clean($chank));
	    $expression_flags {$words[0]} = $words[7];
	    undef @words;
	    }		
    }
    close (EXPRESSION_FLAGS);
}

#------------------------------------------------------------------------------
#read file with occupanicies

my $filesize = -s $in_file; #determine file size in bytes
my $size_counter_step=int($filesize/100);
$filesize = int($filesize/1048576); # filesize in megabytes

print STDERR "\nReading occ column from $in_file file of $filesize MBs. Please wait...\n";

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

my $regex_split_tab='.*\t(.*)';
my $regex_split_newline='\n';
my $regexp_pattern1 = '(..)\t(.*)\t\[(.*)\]';
my $regexp_pattern2 = '^(VIII|VII|VI|V|IV|IX|XIII|XII|XIV|XI|XIX|XVIII|XVII|XVI|XV|X|XX|XXIII|XXII|XXI|III|II|I)\t(\d*)\t(\d*)$';
my $regexp_pattern3 = '(^chr\S{1,2})\s(\d*)\s(\d*)$';
my $regexp_pattern4 = '(^\d*)\s(\d*)$';
my $processed_memory_size = 0;
my $offset=0;

my %occupancy;
my %coords;

while ((my $n = read(INPUT, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <INPUT>;
    }
    my @lines = split(/$regex_split_newline/o, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
	if ($line =~ /$regexp_pattern1/) {
	    $line =~ /$regexp_pattern1/;
	    my $chrom = $1;
	    my $pos = $2; $pos+=0;
	    my $rest = $3;
	    my @rest = split(/, \'\\t\', /, $rest);
	    # normreads
	    #my $plus_val = $rest[3];
	    #my $minus_val = $rest[4];
	    
	    #rawreads
	    my $plus_val = $rest[0];
	    my $minus_val = $rest[1];
	    
	    $chrom =~ s/^0//;
	    if ($invert_strand eq "no") {
		$occupancy{$chrom}{$pos}{'plus'} =  $plus_val;
		$occupancy{$chrom}{$pos}{'minus'} = $minus_val;
	    }
	    else {
		$occupancy{$chrom}{$pos}{'minus'} =  $plus_val;
		$occupancy{$chrom}{$pos}{'plus'} = $minus_val;
	    }
	}
	#elsif ($line =~ /$regexp_pattern2/) {
	#    my $chrom = arabic($1) if isroman($1);
	#    my $pos = $2; $pos+=0;
	#    my $occup_val = $3; $occup_val+=0;
	#    $occupancy{$chrom}{$pos}=$occup_val;
	#    
	#    $input_occ = "yes";
	#}
	elsif ($line =~ /$regexp_pattern3/) {
	    if (!$3) { next; }
	    my $chrom=$1;
	    my $pos = $2; $pos+=0;
	    my $occup_val = $3; $occup_val+=0;
	    $occupancy{$chrom}{$pos}=$occup_val;
	    
	    $input_occ = "yes";
	}
	elsif ($line =~ /$regexp_pattern4/) {
	    my $chrom=$Chromosome;
	    my $pos = $1; $pos+=0;
	    my $occup_val = $2; $occup_val+=0;
	    $occupancy{$chrom}{$pos}=$occup_val;
	    
	    $input_occ = "yes";
	}

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



my @average_occ_freq_distr;
my @output_array;
my @output_array2;
my @non_zero_counter = (0) x ($delta_1+$delta_2+1);
print STDERR "Calculate occupancy\n";
my $genes_counter=0;
my @splice_array;

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

my @length_array;

foreach my $chrom (sort { $a<=>$b || $a cmp $b } keys %occupancy) {
    my @occ_array_plus;
    my @occ_array_minus;
    
    for my $pos  ( sort { $a<=>$b || $a cmp $b } keys %{ $occupancy{$chrom} } ) {
	if ($input_occ eq "no") {
	    $occ_array_plus[$pos] = $occupancy{$chrom}{$pos}{"plus"};
	    $occ_array_minus[$pos] = $occupancy{$chrom}{$pos}{"minus"};
	}
	else {
	    $occ_array_plus[$pos] = $occupancy{$chrom}{$pos};
	    $occ_array_minus[$pos] = $occupancy{$chrom}{$pos};
	}
    }
    my $chr_size= keys %{ $occupancy{$chrom} };
    print STDERR "\nChr $chrom of $chr_size bases...";
    
    my $work_progress_step = int($#TS_positions/10);
    my $current_progress = $work_progress_step;

    for (my $j=0; $j<= $#TS_positions; $j++) {
	my ($position_start, $position_end, $gene_id, $max_occup_counts, $transcript_length);
	undef @splice_array;
	#increment counter to display work progress
	if($current_progress == $j) {print STDERR ".";$current_progress+=$work_progress_step;}

	# check chromosome Nr.
	$chromosomes[$j] =~ s/E-coli/1/;
	$chromosomes[$j] =~ s/Mt/mi/;
	if ($chrom ne $chromosomes[$j]) { next; }
	if (!$chromosomes[$j]) { next; }
	
	# methylation filter
	if ($apply_methylation_filter eq "yes") {
	    if(($Methylation[$j] <= $methylation_range_left) || ($Methylation[$j] >= $methylation_range_right)) { $removed_by_methylation_filter++; next; }
	}	
	# check for repetetive transcript IDs.
	$gene_id = $GeneIDs[$j];
	if($old_id eq $gene_id) { $removed_byID++; next; }
	else { $old_id=$gene_id; }    

	# remove genes by expression flags
	if ((keys %expression_flags) && ($expression_flags {$gene_id} eq "excluded")) { $removed_byExprFlag++; next; }
     
	 # remove minus-strand
	if ($remove_minus_strand eq "yes") {	
	    if ($strands[$j] eq "minus") { $removed_minus++;  next; }
	    else { $position_start=$TS_positions[$j]; $position_end=$TE_positions[$j];}
	}
	
	if (!$strands[$j]) {
	    #code
	    $strands[$j]="plus";
	}
	
	if($old_strand eq $strands[$j]) { $ignore_overlap="no"; }
	else { $ignore_overlap="yes"; $old_strand=$strands[$j]; }

	
	my $check_TSS_position;
	# get transcription start position for plus and minus strand
	if ($strands[$j] eq "minus") {
	    $position_start=$TE_positions[$j];
	    $position_end=$TS_positions[$j];
	    $check_TSS_position=$old_transcript_start-$overlap;
	}
	else {
	    $position_start=$TS_positions[$j];
	    $position_end=$TE_positions[$j];
	    $check_TSS_position=$old_transcript_start+$overlap;
	}
	
	if ($overlap==0) { $ignore_overlap="yes"; }
	    
	# remove string if new TSS is to close to old one and transcript are on different strands
	if ($ignore_overlap eq "no") {
	    if(( $check_TSS_position > $position_start ) && ($strands[$j] eq "plus")) { $removed_overlap++; $ignore_overlap="yes"; next; }
	    elsif(( $check_TSS_position < $position_start ) && ($strands[$j] eq "minus")) { $removed_overlap++; $ignore_overlap="yes"; next; }
	}
	else { $old_transcript_start=$position_start; $ignore_overlap="yes";}
	$transcript_length=abs($position_start-$position_end);    
	
	#initialize array of 0 of 2*$delta+1 size 
	# read part of the array @occ_array; leave 0 if empty
	my ($start_of_region_occ, $end_of_region_occ);
	my ($start_of_region_splice, $end_of_region_splice);
	my ($start_shift,$end_shift);
	my ($old_central_point,$new_central_point);
	
	if ($use_centre eq "yes") {
	    #code
	    $new_central_point=int(($TE_positions[$j]+$TS_positions[$j])/2);
	
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

	    for (my $i=0; $i<= $end_of_region_splice; $i++) {
		my $index=$start_of_region_occ+$i; my $occ_value = $occ_array_plus[$index];
		if($occ_array_plus[$index]) {
		    push(@splice_array, $occ_value);
		    }
		else { push(@splice_array,0); }
	    }
	} else {
	    if ($Cut_tail eq "yes" ) {
		if($strands[$j] eq "plus") {
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
		    
		    for (my $i=0; $i<= $end_of_region_splice+1; $i++) {
			if($occ_array_plus[$start_of_region_occ+$i]) { push(@splice_array, $occ_array_plus[$start_of_region_occ+$i]); }
			else { push(@splice_array,0); }
		    }
		    
		    if ($invert_strand eq "no") {
			if ($end_of_region_splice < $delta_1+$delta_2+1) {
			    #code
			    push(@splice_array,0) for($end_of_region_splice+1..$delta_1+$delta_2+1);
			}
		    }
		    else {		    
			if ($end_of_region_splice < $delta_1+$delta_2+1) {
			    #code
			    push(@splice_array,0) for($end_of_region_splice+1..$delta_1+$delta_2+1);
			}
			my @temp = @splice_array;
			@splice_array = reverse(@temp);
		    }
		    
		}
		elsif ($strands[$j] eq "minus") {
		    # position START and END swapted!!
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
				    
		    $start_of_region_splice=0;
		    
		    for (my $i=0; $i<=$end_of_region_splice+1; $i++) {
			if($occ_array_minus[$start_of_region_occ+$i]) { push(@splice_array, $occ_array_minus[$start_of_region_occ+$i]); }
			else { push(@splice_array,0); }
		    }
				    
		    if ($invert_strand eq "no") {
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
			}
		    }
		    
    
		}
		if ($verbose==1) {
		    print STDERR join("\t","$j: $strands[$j]", $gene_id, $chrom , $chromosomes[$j], "TSS: $position_start", "TTS: $position_end", "TSS-delta1: $start_of_region_occ","TTS/TSS+delta2: $end_of_region_occ","$start_of_region_splice","$end_of_region_splice"),"\t";
		    print STDERR "$#splice_array\n";
		}
	    }
	    
	    elsif ($Cut_tail eq "no" ) {
		$start_of_region_splice=0;
		$end_of_region_splice=$delta_1+$delta_2+1;
    
		if ($strands[$j] eq "plus") {
		    
		    $start_of_region_occ = $position_start-$delta_1;
		    $end_of_region_occ = $position_start+$delta_2+1;
		    
		    for (my $i=0; $i<= $end_of_region_splice; $i++) {
			if($occ_array_plus[$start_of_region_occ+$i]) { push(@splice_array, $occ_array_plus[$start_of_region_occ+$i]); }
			else { push(@splice_array,0); }
		    }
		    if ($verbose==1) {
			print STDERR join("\t","$j: $strands[$j]", $gene_id,  $chrom , $chromosomes[$j], "TTS: $position_start", "TSS: $position_end", "TTS/TSS-delta2: $start_of_region_occ","TTS+delta1: $end_of_region_occ","$start_of_region_splice","$end_of_region_splice"),"\t";
		    }
		}
		elsif ($strands[$j] eq "minus") {
		    $start_of_region_occ = $position_start-$delta_2;
		    $end_of_region_occ= $position_end+$delta_1+1;
		    
		    for (my $i=0; $i<=$end_of_region_splice; $i++) {
			if($occ_array_minus[$start_of_region_occ+$i]) { push(@splice_array, $occ_array_minus[$start_of_region_occ+$i]); }
			else { push(@splice_array,0); }
		    }
		    
		    my @temp = @splice_array;
		    @splice_array = reverse(@temp);
		    
		    if ($verbose==1) {
			print STDERR join("\t","$j: $strands[$j]", $gene_id, "\t",  $chrom , $chromosomes[$j], "TTS: $position_start", "TSS: $position_end", "TTS/TSS-delta2: $start_of_region_occ","TTS+delta1: $end_of_region_occ","$start_of_region_splice","$end_of_region_splice"),"\t";
		    }
		} 
	    }

	}
	
	# check for artefacts
	$max_occup_counts = max(@splice_array);
	if (($max_occup_counts > $upper_occup_threshold) or ($max_occup_counts < $lower_occup_threshold)){ $removed_threshold++; next; }
	 	
	# normalize
	if ($GeneLengthNorm eq "yes") {
	    #code
	    my $total_reads_per_transcript;
	    $total_reads_per_transcript += $_ for @splice_array;
	    my $norm_factor=$total_reads_per_transcript/$transcript_length;
	    if ((!$total_reads_per_transcript) or ($total_reads_per_transcript == 0)) { @splice_array = (0) x ($delta_1+$delta_2+1); }
	    else {
		my @temp_array = map { $_ / $norm_factor } @splice_array;
		undef @splice_array;
		@splice_array=@temp_array;
	    }
	}
	
	
	if($library_size_normalization eq "yes") {
	    my $norm_factor=$library_size;
	    if ((!$norm_factor) or ($norm_factor == 0)) { @splice_array = (0) x ($delta_1+$delta_2+1); }
	    else {
		my @temp_array = map { $_ / $norm_factor } @splice_array;
		undef @splice_array;
		@splice_array=@temp_array;
	    #print STDERR $norm_factor,"\t";
	    }
        }

	
	if($apply_DivSum_normalization eq "yes") {
	    my $norm_factor += $_ for @splice_array;
	    if ($norm_factor == 0) { @splice_array = (0) x ($delta_1+$delta_2+1); }
	    else {
		my @temp_array = map { $_ / $norm_factor } @splice_array;
	    undef @splice_array;
	    @splice_array=@temp_array;
	    #print STDERR $norm_factor,"\t";
	    }
        }
	
	my $temp_string;
	if ($calc_score eq "yes") {
	    my $RPM = sum(@splice_array)/(1000000*$library_size_normalization);
	    $temp_string = join("\t", $gene_id, $transcript_length, $RPM);
	}
	else {
	    $temp_string = join("\t", $gene_id, $transcript_length, @splice_array);
	}

	push (@output_array, $temp_string);
	
	$temp_string = join("\t", $gene_id, $chromosomes[$j], $strands[$j], $start_of_region_occ, $end_of_region_occ, $end_of_region_occ-$start_of_region_occ , $TS_positions[$j], $TE_positions[$j], $transcript_length, $#splice_array);
	push (@output_array2, $temp_string);
	
	#increment each element of positioning array
	for (my $i=0; $i<= $#splice_array; $i++) {
	    $average_occ_freq_distr[$i] += $splice_array[$i];
	    $non_zero_counter[$i]++;
	}
	$genes_counter++;
	push(@length_array, $transcript_length);
	#print STDERR "$temp_string\n";
    }
    
    undef @occ_array_plus;
    undef @occ_array_minus;

}

my @results;
my @TrasncriptsNr_perBase;
for (my $i=0; $i<=$#average_occ_freq_distr; $i++) {
    my $norm_factor=checkarray(\@length_array,$i);
    push(@TrasncriptsNr_perBase,$norm_factor);
}

# generate average occupancy distribution
if ( $normalize eq "yes" ) {
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
if ( $PerBaseNorm eq "yes" ) {
	for (my $k=0;$k<=$#average_occ_freq_distr; $k++ ) {
	    if ($verbose==1) {
		print STDERR "$k:\t$average_occ_freq_distr[$k]\t$TrasncriptsNr_perBase[$k]\t";
	    }
	    if ($TrasncriptsNr_perBase[$k] == 0) {
		#code
		$TrasncriptsNr_perBase[$k]=1;
	    }
	    
	    my $normalized_occup_summ = $average_occ_freq_distr[$k] / $TrasncriptsNr_perBase[$k];
	    if ($verbose==1) {
		print STDERR "$normalized_occup_summ\n";
	    }
	    push(@results,$normalized_occup_summ);
	}
}
else { @results = @average_occ_freq_distr; }


# writing outputs
if ($dont_save_aligned eq "no") {
    @output_array = grep /\S/, @output_array;
    open (OCCUP_File, ">$out_path1") or die "can't open file $out_path1 for writting: $!";
    print OCCUP_File join("\n", @output_array),"\n";
    close (OCCUP_File);
    
    $out_path1 =~ s/\.txt/\.tab/;
    open (OCCUP_File, ">$out_path1") or die "can't open file $out_path1 for writting: $!";
    print OCCUP_File join("\t", "Gene", "Chromosomes", "strand", "start_of_region_occ", "end_of_region_occ", "region length" , "TSS", "TTS", "transcript_length", "array length"),"\n";
    print OCCUP_File join("\n", @output_array2),"\n";
    close (OCCUP_File);

}

open (AveragedFile, ">$out_path2") or die "can't open file $out_path2 for writting: $!";
print AveragedFile join("\n", @results),"\n";
close(AveragedFile);

print STDERR "done\n$#output_array transcription starts has been identified.\n\n",
"Removed transcripts:\n",
"minus strand: $removed_minus\n",
"TSS overlap: $removed_overlap\n",
"threshold: $removed_threshold\n",
"methylation: $removed_by_methylation_filter\n",
"identical (by ID): $removed_byID\n",
"removed by expression fllag (excluded after CW): $removed_byExprFlag\n",
"$removed_by_chromosome from $#TS_positions transcripts removed (other chromosomes)\n\n";

open(LOG, ">>$run_log_path") or die "can't uppend to a log file: $run_log_path for writting: $!";
print LOG "list\t$genes_counter\n";
close (LOG);


undef @LIST_array;
undef @TS_positions;
undef @TE_positions;
undef @strands;
undef @chromosomes;
undef @GeneIDs;
undef @Expression;
undef @splice_array;
undef @non_zero_counter;
undef @Methylation_col1;
undef @Methylation_col2;
undef @Methylation;


$tm = localtime;
my $stop_sec = $tm -> [0];
my $stop_min = $tm ->[1];
my $stop_hour = $tm ->[2];;
my $stop_time = time();
my $message;
$duration = $stop_time-$start_time;

$message = "\nStarted:\t$start_hour:$start_min:$start_sec\nnow:\t$stop_hour:$stop_min:$stop_sec\nduration:\t$duration sec.\n";

print STDERR "$message\nJob finished!\nBye!\n\n";

exit;



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
    my ($column_number, $column_name, $array_ref) = @_;
    my @array = @{$array_ref};

    my (@column, @string);
    # read column of interest to the memory
    for(my $j=0; $j <= $#array ; $j++ )
     {
	    push (@string, split("[\t\r\f\n\,]",$array[$j]));
	    if ($column_number > $#string) {push (@column,undef);}
	    else {push (@column, $string[$column_number]);}
	    undef @string;
     }
    print STDERR join("\n", $column_name, "_________", $column[0], $column[1], $column[2]), "\n";
    return (@column);
}

#-------------------------------------------------------------
sub clean {

my $text = shift;

$text =~ s/\r//g;
$text =~ s/\n//g;
return $text;
}

sub max {
  my $max = $_[0];
  for ( @_[ 1..$#_ ] ) {
    $_ = 0 if !$_;
    $max = $_ if $_ > $max; }
  return($max);
}

sub min {
  my $min = $_[0];
  for ( @_[ 1..$#_ ] ) {
    $_ = 0 if !$_;
    $min = $_ if $_ < $min; }
  return($min);
}

# catch exceptions
sub try(&) { eval {$_[0]->()} };
sub catch(&) { $_[0]->($@) if $@ }


__DATA__
###==================================================================================================
### Calculates aggregate profile of sequencing read density around genomic regions
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### aggregate_profile.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 19 July 2016
###==================================================================================================


 Note: Program runs in command-line mode only
 ____________________
 Usage instruction:
 ____________________
perl -w aggregate_profile.pl -input=OccupancyFileName -regions=GenomicRegionsFileName -expression=GeneExpressionFileName -aligned=[N] -average_aligned=[N] -path2log= -region_start_column=[N] -region_end_column=[N] -strand_column=[N] -chromosome_col=[N] -GeneId_column=[N] -Expression_columnID=[N] -Methylation_columnID=[N] -Methylation_columnID2=[N] -upstream_delta=[N] -downstream_delta=[N] -upper_threshold=[N] -lower_threshold=[N] -Methylation_threshold=[N|N-N] -overlap=[N] -library_size=[N] -remove_minus_strand | -ignore_strand | -fixed_strand=[plus|minus] -invert_strand -input_occ -score -dont_save_aligned -Cut_tail -chromosome=chrN -AgregateProfile -GeneLengthNorm -LibsizeNorm -PerBaseNorm -useCentre
  -use_default --verbose --help. The program’s help mode contains detailed explanations.

Parameters explained:
 _______________________
-input="file.occ": file with the occupancy data for a given chromosome;
-regions="regins_annotation.txt": path to a regions annotation file;
-expression="GeneExpression.txt": tab-delimited text file, with RPKM values;
-GeneId_column=[N]: number of column with gene ID (default column Nr. 0);
-chromosome_col=[N]: column with chromosome Nr. (default column Nr. 6);
-strand_column=[N]: number of columns with strand (default column Nr. 7);
-region_start_column=[N]: default column Nr. 8;
-region_end_column=[N]: default column Nr. 9;
-Expression_columnID=[N]: default column Nr. 14;
-Methylation_columnID=[N]: default column Nr. 2;
-Methylation_columnID2=[N]: default column Nr. 3;
-upstream_delta=[N]: number of base pairs upstream from aligned starts considered in calculations. Default value is 100;
-downstream_delta=[N]: number of base pairs downstream from aligned starts considered in calculations. Default value is 1500;
-upper_threshold=[N]: set upper occupancy threshold to remove mapping/sequencing artefacts. Default value is 10000;
-lower_threshold=[N]: set lower occupancy threshold to avoid ORFs with low coverage. Default value is 0;
-Methylation_threshold=[N]: set the threshold for the methylation value. Keep only sites with methylation within the range. Default range: 50-100%.
-Methylation_threshold=[N1-N2];
-overlap=[N]: when TSS for differnt transcript IDs occurs within a range of [overlap] thay considered as duplicates and removed from calcualtions. Deafult value is 100;
-library_size=[N]: sequencing library size (to use with -LibsizeNorm to normalize to a sequencing library size);
-dont_save_aligned: do not save the generated aligned matrix;
-score: convert occupancy to RPM values and use it to generate aggregate profile;
 -chromosome=chrN: limit analysis to regions derived from specified chromosome only;
 -useCentre: Use middle of the region for alignment instead of start of the region;
 
 Normalization settings:
 _______________________
 -AgregateProfile: generate average occupancy profile (one column text table);
  -GeneLengthNorm: normalize aggregate profile of each region to it length (generate ametagene);
  -LibsizeNorm: normalize to a sequencing library size;
  -PerBaseNorm: perform per-base normalization of aggregate profile to compensate the transcripts length difference;

 Strand-related settings:
 ________________________
 -remove_minus_strand: do not use regions on minus strand for alignment and aggregate profile calculations;
 -ignore_strand: do not use strand information to do alignment;
 -fixed_strand=[plus|minus]: input reads belongs only to either plus or minus strand;
 -invert_strand: invert strand from annotation to opposite (use for aligning to the end of genomic regions);
 
 OUTPUT files:
 ____________________
 -aligned="aligned.txt": path to the output file with aligned to a region start
 -average_aligned="aligned_average.txt": path to the output file with average occupancy profile

 LOG file:
 ____________________
 -path2log="run.log": path to a tab-delimited text file, containing number of analyzed transcripts per chromosome/list 
  
 Additional parameters:
 ____________________
 -use_default: use default parameters and paths to the files.
 --help: display this message
 
==================================================================================================
