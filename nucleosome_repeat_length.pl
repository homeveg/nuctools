#!/usr/local/bin/perl
###
###==================================================================================================
### Calculates frequency of nucleosome-nucleosome distances to determine the nucleosome repeat length
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### nucleosome_repeat_length.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================


use strict "vars";
use Config;
use Time::localtime;
use Time::Local;


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

$in_file = "chr1.bed";
$out_path1 = "nuc-nuc_chr1.txt";
my $apply_filter_flag=0;
my $piles_filtering_threshold=20;
my $fix_pile_size = "no";

# perl -w new_nuc-nuc_distance_filter.pl -input="chr9.bed" -output="nuc-nuc_ch9_filtered.txt" -delta=3000 -filtering_threshold=20 -apply_filter

#read arguments from command line
if (@ARGV != 0) {
    foreach my $comand_line_flag (@ARGV) {
	if ($comand_line_flag =~ /-input=(.*)/i) { $in_file = $1; }
	if ($comand_line_flag =~ /-output=(.*)/i) { $out_path1 = $1; }
        if ($comand_line_flag =~ /-delta=(.*)/i) { $delta = $1; }
	if ($comand_line_flag =~ /-pile=(.*)/i) { $pile = $1; }
	if ($comand_line_flag =~ /-fix_pile_size/i) { $fix_pile_size = "yes"; }
        if ($comand_line_flag =~ /-use_default/i) { print STDERR "using default values: \n"; }
        if ($comand_line_flag =~ /-apply_filter/i) { $apply_filter_flag=1; }
        if ($comand_line_flag =~ /-filtering_threshold=(.*)/i) { $piles_filtering_threshold = $1; }
	
	
	#if ($comand_line_flag =~ /-help/i) {  print STDERR "help\n";	}
    }
}
else { warn "please provide command line options!\n:",
      "perl -w nuc-nuc_distance.pl -input -output -delta "; exit;}

# Display input parametrs
print STDERR "======================================\n";
print STDERR "Started:\t$start_hour:$start_min:$start_sec\n";
print STDERR "======================================\n";
print STDERR "in file:",$in_file, "\n";
print STDERR "out file:",$out_path1, "\n";
print STDERR "delta: ",$delta, "\n";
print STDERR "pile: $pile\n";
print STDERR "select only fix pile size: $fix_pile_size\n";
print STDERR "apply_filter_flag: $apply_filter_flag\n";
print STDERR "filtering threshold: $piles_filtering_threshold\n";
#print STDERR "reads threshold: ",$reads_threshold, "\n";
#print STDERR "using fast files load (nor for files bigger than 2Gb: $use_mmap_perlIO\n";
print STDERR "======================================\n";


my @occ_array=();

#read first file with occupanicies
#read file with by 4kb chanks

print STDERR "Step 1 (",time()-$timer1," sec.): reading $in_file file...\n";

@occ_array=();
my $BUFFER_SIZE = 1024*4;

# open original file
open(INPUT, $in_file) or die "error: $in_file cannot be opened\n";
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
while ((my $n = read(INPUT, $buffer, $BUFFER_SIZE)) !=0) {
    if ($n >= $BUFFER_SIZE) {
    $buffer .= <INPUT>;
    }
    my @lines = split(/$regex_split_newline/o, $buffer);
    # process each line in zone file
    foreach my $line (@lines) {
	chomp($line);
        my @newline1=split(/\t/, $line);
        my $start_nuc=$newline1[1];
        push(@starts, $start_nuc);

        $string_counter++;
	if ($start_nuc>0) {$not_zero_counter++;}
    }
$processed_memory_size += $n;
$offset += $n;
if(int($processed_memory_size/1048576)>= $filesize/10) {
    print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
    }
undef @lines;
$buffer = "";
}

close(INPUT) or die $!; 
print STDERR $filesize, " Mbs processed in ", time()-$timer2, " seconds.\n$not_zero_counter non zero counts, ",$#starts+1," lines\n\n";

# remove empty strings
$timer2= time();
print STDERR "- cleaning from empty strings (if any)...";
@starts = grep /\S/, @starts;
print STDERR "done in ", time()-$timer2, " seconds. ",$#starts+1," strings left\n";

# sort nucleosome positions according to a start_nuc
$timer2= time();
print STDERR "- sorting...";
my @sorted = sort {$a <=> $b} @starts;
print STDERR "done in ", time()-$timer2, " seconds.\n";
$timer2= time();

# remove nucleosomoes without repeat ($pile>1)
my @temp;
if ($pile>1) {
    $timer2= time();
    if($fix_pile_size eq "yes") {
	print STDERR "- select only pile=3...";
    }
    else {
	print STDERR "- removing un-piled nucleosomes...";
    }
    
    my @only_piled;
    my $pile_counter=0;

    for (my $i=1; $i<=$#sorted; $i++) {
	if (!@temp) { push(@temp,$sorted[$i-1]); }
	if ($sorted[$i-1] == $sorted[$i]) {
	    push(@temp,$sorted[$i]);
	    $pile_counter++;
	}
	elsif ($pile_counter < $pile) {
	    undef @temp;
	    $pile_counter=0;
	}
	elsif (($sorted[$i-1] != $sorted[$i]) && ($#temp>0)) {
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
    undef @sorted;
    @sorted = grep /\S/, @only_piled;
    print STDERR "done in ", time()-$timer2, " seconds. ",$#sorted+1," strings left\n";
}

#apply local pile filter: at one position should be less than threshold nucleosom starts
if($fix_pile_size eq "no") {
    $timer2= time();
    undef @temp;
    if ($apply_filter_flag==1) {
	$timer2= time();
	print STDERR "- apply local pile filter: removing reads in the pile above $piles_filtering_threshold ...";
	my @piled_under_threshold;
	my $pile_counter=0;
    
	for (my $i=1; $i<=$#sorted; $i++) {
	    if (!@temp) { push(@temp,$sorted[$i-1]); }
	    if ($sorted[$i-1] == $sorted[$i]) {
		push(@temp,$sorted[$i]);
		$pile_counter++;
	    }
	    elsif ($pile_counter >= $piles_filtering_threshold) {
		push @piled_under_threshold, @temp[0..$piles_filtering_threshold];
		undef @temp;
		$pile_counter=0;
	    }
	    elsif (($sorted[$i-1] != $sorted[$i]) && ($#temp>0)) {
		push @piled_under_threshold, @temp;
		undef @temp;
		$pile_counter=0;
	    }
    
	}
	undef @sorted;
	@sorted = grep /\S/, @piled_under_threshold;
	print STDERR "done in ", time()-$timer2, " seconds. ",$#sorted+1," strings left\n";
    }
}


print STDERR "Step 2: Start analysis...\n";
$timer2= time();    
my @output_array = (0) x ($delta+1); #initialize array of 0 of $delta+1 size
my $first_itteration=0;
my $counter_step=int($#sorted/100);
my $progress_counter=$counter_step;

for (my $i=0; $i<$#sorted; $i++) {
    #read read start
    my $nuc=$sorted[$i];
    # calcualte maximum index shift
    my $max_delta_index=5*$delta;
    # check if incremented index exceeds reads array length and correct it if necessary
    if ($i+$max_delta_index>=$#sorted) { $max_delta_index = $#sorted-$i; }
    for (my $n=1; $n<=$max_delta_index ; $n++) {
	if (!$sorted[$i+$n]) { last; }
        my $nuc_plus=$sorted[$i+$n];
        my $delta_nuc_starts = $nuc_plus-$nuc;
	if (($delta_nuc_starts>$delta) or ($delta_nuc_starts<0)) { last; }
        else { $output_array[$delta_nuc_starts]++; }
    }
    #increment counter to display work progress
    if ( $progress_counter == $i ) {
	if ($first_itteration==0) { my $approx_duration = (time()-$timer2)*100/60; print STDERR "work will be finished in approximately $approx_duration minutes\n"; $first_itteration=1; }
	#print STDERR ".";
	print STDERR join (" ", @output_array[0..5],"..",@output_array[145..150]), "\n";
	$progress_counter+=$counter_step;
    }
}
print STDERR "\ndone\n"; 

@output_array = grep /\S/, @output_array;

print STDERR "- saving resuts to $out_path1...";
open (OCCUP_File, ">$out_path1") or die "can't open file $out_path1 for writting: $!";
print OCCUP_File join("\n", @output_array);
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