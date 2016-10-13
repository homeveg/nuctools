#!/usr/local/bin/perl

###==================================================================================================
### Calculates genome-wide occupancy based on the bed file with sequencing reads
### (c) Yevhen Vainshtein, Vladimir Teif
### 
### bed2occupancy_average.pl
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================

use IO::Dir;

use strict "vars";
my $usage = "$0 -input=<in.bed> -output=<out.occ> -chromosome_col=<column Nr.> -start_col=<column Nr.> -end_col=<column Nr.> -strand_col=<column Nr.> -window=<running window size> -consider_strand -ConvertAllInDir\n";
my ($infile,$outfile);

#default columns

my $flag = 0;
my $ConvertAllInDir=0;

my $start_col=1;
my $end_col=2;
my $strand_col=3;
my $chromosome_col=0;
my $running_window=100;
my $region_start=0;
my $region_end;

if (@ARGV != 0) {
    foreach my $comand_line_flag (@ARGV) {
	if ($comand_line_flag =~ /-input=(.*)/i) { $infile = $1; }
	if ($comand_line_flag =~ /-output=(.*)/i) { $outfile = $1;}
	if ($comand_line_flag =~ /-consider_strand/i) { $flag = 1;}
	if ($comand_line_flag =~ /-ConvertAllInDir/i) { $ConvertAllInDir = 1;}
	if ($comand_line_flag =~ /-start_col=(.*)/i) { $start_col = $1; }
	if ($comand_line_flag =~ /-end_col=(.*)/i) { $end_col = $1; }
	if ($comand_line_flag =~ /-strand_col=(.*)/i) { $strand_col = $1; }
	if ($comand_line_flag =~ /-window=(.*)/i) { $running_window = $1; }
	if ($comand_line_flag =~ /-region_start=(.*)/i) { $region_start = $1; }
	if ($comand_line_flag =~ /-region_end=(.*)/i) { $region_end = $1; }
	
	if ($comand_line_flag =~ /-chromosome_col=(.*)/i) { $chromosome_col = $1; }
    }
}
else { die $usage; }

# Display input parametrs
print STDERR "======================================\n";
print STDERR "input file or folder: ",$infile, "\n";
print STDERR "output file (if applicable): ", $outfile, "\n";
print STDERR "================ flags ==============\n";
print STDERR "convert all files in the dirrectory: ",$ConvertAllInDir, "\n";
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



if ($ConvertAllInDir == 0) {     Bed2Occup($infile, $outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag); } 
elsif ($ConvertAllInDir == 1) {
    # process each *.bed file in the folder
    my (%dir, @dir_list, @text_list);
    my $start_dir = $infile;
    tie %dir, IO::Dir, $start_dir;
    foreach (keys %dir) { push (@dir_list, $_); }
    my $default_name = "bed_files.lst";

    foreach my $file (@dir_list) {
        if ($file =~ m/.*\.bed$/) {
            $outfile = $file;
            $outfile =~ s/(.*)\.bed$/$1\.occ/;
            Bed2Occup($file, $outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag);           
        }
    }
}

exit;

#--------------------------

sub Bed2Occup {
    my ($infile, $outfile, $chromosome_col, $start_col, $end_col, $strand_col, $flag) = @_;
    
    open(IN, $infile) || die "Can't open $infile for reading!\n";
    open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";
    
    my $buffer = "";
    my $sz_buffer = 0;
    my $timer2 = time();
    # counter for the markers we see
    my $marker_count = 0;
    
    my $regex_split_newline='\n';
    
    my $filesize_in_bytes = -s $infile; #determine file size in bytes
    my $size_counter_step=int($filesize_in_bytes/10);
    my $filesize = int($filesize_in_bytes/1048576); # filesize in megabytes
    
    print STDERR "Reading $infile file of $filesize MBs. Please wait...\n";
    my $processed_memory_size = 0;
    my $offset=0;
    my $not_zero_counter=0;
    my $string_counter=0;
    my $BUFFER_SIZE = 1024;
    my @occup=();
    my $old_coordinate=1;
        
    while ((my $n = read(IN, $buffer, $BUFFER_SIZE)) !=0) {
        if (($n >= $BUFFER_SIZE) or (($n == $filesize_in_bytes))) {
            $buffer .= <IN>;
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
            print STDERR int($offset/1048576), " Mbs processed in ", time()-$timer2, " seconds.\n"; $processed_memory_size=0;
            #last;
            }
        undef @lines;
        $buffer = "";
    }
    close(IN);
    
    print STDERR "saving normalized occupancy data.\nThere are ", $#occup+1, " lines in new output file\n";
    open(OUT, ">$outfile") || die "Can't open $outfile for writing: $!\n";
    $timer2 = time();
 
    #initialize running average for first 2*$running_window+1 nucleotides
    my ($sum, $normalized_occupancy);
    print STDERR "calculating and printing normalized occupancy with a running window +/-",$running_window,"...";
    
    # modify running average by shifting
    if (!$region_end) {	$region_end=$#occup; }
    for (my $i=$region_start; $i<$region_end; $i+=$running_window) {
	$sum=0;
	for (my $j=$i; $j<=$i+$running_window; $j++) {
	    $sum+=$occup[$i];
	}
	my $average = $sum/$running_window;
	
	if ($average==0) {
	    print OUT join ("\t", $i+$running_window,0),"\n";
	}
	else {
	    $normalized_occupancy = $occup[$i+int($running_window/2)]/$average;
	    print OUT join ("\t", $i+$running_window,$average),"\n";	    
	}
	undef $sum;
    }
    print STDERR "done\n";

    close(OUT);    
    print STDERR "Job finished!\n$outfile generated in ", time()-$timer2, " seconds.\nGood bye!\n";

}
