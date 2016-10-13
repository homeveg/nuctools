#!/usr/bin/perl

###==================================================================================================
### Merge 2 tables by the columns with unique IDs
### (c) Yevhen Vainshtein
### 
### merge2tabs.pl
###
### NucTools 1.0
###==================================================================================================
###
### last changed: 17 July 2016
###==================================================================================================

use strict;

my $usage = "$0 -table1=\"path to table 1\" -table2=\"path to table 2\" -output=\"path to output\" -colID_tab1=<Column number from table 1> -colID_tab2=<Column number from table 1> -header\n";

my ($path_tab1, $path_tab2, $path2output, $colID_tab1,$colID_tab2);
#read arguments from command line
if (@ARGV != 0) {
    foreach my $comand_line_flag (@ARGV) {
	#input files
	if ($comand_line_flag =~ /-table1=(.*)/i) { $path_tab1 = $1; } 
        if ($comand_line_flag =~ /-table2=(.*)/i) { $path_tab2 = $1; } 
        if ($comand_line_flag =~ /-output=(.*)/i) { $path2output = $1; } 
        if ($comand_line_flag =~ /-colID_tab1=(\d*)/i) { $colID_tab1 = $1; } 
        if ($comand_line_flag =~ /-colID_tab2=(\d*)/i) { $colID_tab2 = $1; } 
    }
} else {
    warn $usage;
    exit;
}

print STDERR "table 1:",$path_tab1, "\n";
print STDERR "table 2:",$path_tab2, "\n";
print STDERR "column with ID, table 1:",$colID_tab1, "\n";
print STDERR "column with ID, table 2:", $colID_tab2, "\n";
print STDERR "save merged table to a file: ",$path2output, "\n";

my (@array_tab1,@array_tab2);
open(INPUT, $path_tab1) or die "$path_tab1 cannot be opened: $!\n";
while (<INPUT>) { for my $chank  (split/\r\n/) { my $text = clean($chank); push(@array_tab1, $text); } }
close(INPUT);

open(INPUT, $path_tab2) or die "$path_tab2 cannot be opened: $!\n";
while (<INPUT>) { for my $chank  (split/\r\n/) { my $text = clean($chank); push(@array_tab2, $text); } }
close(INPUT);

my ($array_tab1_ref, $array_tab2_ref,$match_array1_ref,$match_array2_ref) = compare_2_arrays($colID_tab1, $colID_tab2, \@array_tab1, \@array_tab2);

my @joined_arrays = join_arrays($match_array1_ref,$match_array2_ref);

open(OUTPUT, ">$path2output") or die "$path2output cannot be opened for writing: $!\n";

print OUTPUT join("\n", @joined_arrays), "\n";
close(OUTPUT);
print STDERR "\nfinished\n";
exit;

#-------------------------------------------------------------
#-------------- match 2 arrays--------------------------------
sub compare_2_arrays {
    my ($selected_column1, $selected_column2, $array_ref1, $array_ref2) = @_;
    my @array1 = @{$array_ref1};
    my @array2 = @{$array_ref2};
    
    my $nr_of_cols1 = () = $array1[0] =~ /\t/g;
    my $nr_of_cols2 = () = $array2[0] =~ /\t/g;
    my $empty_string_array1 = "\t"x$nr_of_cols1;
    my $empty_string_array2 = "\t"x$nr_of_cols2;
    
    my $counter1=1;
    my $counter2;
    
    my (@match_array1, @match_array2);

    my @sel_column_array1 = Read_column($selected_column1, $array_ref1);   
    my @sel_column_array2 = Read_column($selected_column2, $array_ref2);   
        
    my $progress_step = int($#array1/20); my $progress_counter=0;

    while ($counter1<=$#array1) {
	$counter2=1;
	my $ID_A = $sel_column_array1[$counter1]; #print STDERR join("\t",$counter1,$ID_A);
	while ($counter2<=$#array2) {
	    my $ID_B = $sel_column_array2[$counter2];
            #print STDERR "\t\t",join("\t",$counter2,$ID_B),"\n";
	    if ("$ID_A" eq "$ID_B") {
		push(@match_array1, $array1[$counter1]);
		push(@match_array2, $array2[$counter2]);
		splice (@array1, $counter1, 1, $empty_string_array1);
		splice (@array2, $counter2, 1, $empty_string_array2);
		last;
	    }
	    $counter2 += 1;
	}
        #print STDERR "\n";
	$counter1 += 1;    
	$progress_counter++;
	if($progress_counter==$progress_step) {print STDERR "."; $progress_counter=0;}
    }
    my (@rest_array1,@rest_array2);
    foreach my $element (@array1) {
	if ($element) {push(@rest_array1,$element);}
    }
    foreach my $element (@array2) {
	if ($element) {push(@rest_array2,$element);}
    }
    
    return (\@array1, \@array2,\@match_array1,\@match_array2,\@rest_array1,\@rest_array2);
}

#-------------------------------------------------------------
#------------------ join 2 arrays of the same size -----------
sub join_arrays {
    my ($array_ref1, $array_ref2) = @_;
    my @array1 = @{$array_ref1};
    my @array2 = @{$array_ref2};
    my @array;
    for (my $i=0; $i <= $#array1; $i++) {
	push (@array, "$array1[$i]\t\t$array2[$i]");
    }
    return(@array);
}


#-------------------------------------------------------------
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
    return (@column);
}

#-------------------------------------------------------------
#------------------ Read specified string --------------------
sub Read_string
{
    my ($string_number, $array_ref) = @_;
    my @array = @{$array_ref};

push (my @string, split ("\t",$array[$string_number]));
return (@string);
}

#-------------------------------------------------------------
#------------------ clean line endings -----------------------
sub clean {

my $text = shift;

$text =~ s/\r//g;
$text =~ s/\n//g;
return $text;
}
