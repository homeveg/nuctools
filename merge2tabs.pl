#!/usr/bin/perl

=head1 NAME

merge2tabs.pl - Merge 2 tables by the columns with unique IDs

=head1 SYNOPSIS

perl -w merge2tabs.pl --table1=<path to table 1> --table2=<path to table 2> --output=<path to merged> [--colID_tab1=<column Nr. tab1> --colID_tab2=<column Nr. tab2> --gzip --help]

 Required arguments:
    --table1 | -t1       path to table 1
    --table2 | -t2       path to table 2
    --output | -out      path to output merged table

 Options:
    --colID_tab1 | -c1   column Nr. with unique IDs, table 1 (default: 0)
    --colID_tab2 | -c2   column Nr. with unique IDs, table 2 (default: 0)
	
	--gzip               compress resulting merged table with gzip (when found in the system path)
	--help | h           Help
	
 Example usage:
 
    perl -w merge2tabs.pl --table1=table1.txt.gz --table2=table2.txt.gz --output=merged.txt.gz --colID_tab1=1 --colID_tab2=2
	
	OR
	
    perl -w merge2tabs.pl -t1 table1.txt.gz -t2 table2.txt.gz -out merged.txt.gz -c1=1 -c2=2
    
=head1 DESCRIPTION
 
=head2 NucTools 1.0 package.

 NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data

=head2 average_replicates.pl

 extend_PE_reads.pl Merge 2 tables by the columns with unique IDs

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

use strict;
use Getopt::Long;
use Pod::Usage;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

my ($path_tab1, $path_tab2, $path2output);
my $colID_tab1 = 0;
my $colID_tab2 = 0;

my $needsHelp;
my $useGzip;

my $options_okay = &Getopt::Long::GetOptions(
	'table1|t1=s' => \$path_tab1,
	'table2|t2=s'   => \$path_tab2,
	'output|out=s' => \$path2output,
	'colID_tab1|c1=s'   => \$colID_tab1,
	'colID_tab2|c2=s'   => \$colID_tab2,

	'gzip'      => \$useGzip,
	'help|h'      => \$needsHelp
);

# Check to make sure options are specified correctly and files exist
&check_opts();

# check if GZIP is loaded
if ( ((!$ModuleGzipIsLoaded) or (!$ModuleGunzipIsLoaded)) and ($useGzip) ) {
	print STDERR "Can't work with GZIP: IO::Compress::Gzip is not on PATH\n";
	exit;
}
elsif ( (($ModuleGzipIsLoaded) and ($ModuleGunzipIsLoaded)) and ($useGzip) ) {
	print STDERR "ZGIP support enabled\n";
}
else {
	print STDERR "ZGIP support disabled\n";
}
print STDERR "table 1:",$path_tab1, "\n";
print STDERR "table 2:",$path_tab2, "\n";
print STDERR "column with ID, table 1:",$colID_tab1, "\n";
print STDERR "column with ID, table 2:", $colID_tab2, "\n";
print STDERR "save merged table to a file: ",$path2output, "\n";

my (@array_tab1,@array_tab2);
# open occupancy file

my $inFH1;
if ( $path_tab1 =~ (/.*\.gz$/) ) {
	$inFH1 = IO::Uncompress::Gunzip->new( $path_tab1 )
	or die "IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip::GunzipError\n";
}
else { open( $inFH1, "<", $path_tab1 ) or die "error: $path_tab1 cannot be opened:$!"; }
while (<$inFH1>) { for my $chank  (split/\r\n/) { my $text = clean($chank); push(@array_tab1, $text); } }
close($inFH1);


my $inFH2;
if ( $path_tab2 =~ (/.*\.gz$/) ) {
	$inFH2 = IO::Uncompress::Gunzip->new( $path_tab2 )
	or die "IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip failed::GunzipError\n";
}
else { open( $inFH1, "<", $path_tab2 ) or die "error: $path_tab2 cannot be opened:$!"; }
while (<$inFH2>) { for my $chank  (split/\r\n/) { my $text = clean($chank); push(@array_tab2, $text); } }
close($inFH2);

my ($array_tab1_ref, $array_tab2_ref,$match_array1_ref,$match_array2_ref) = compare_2_arrays($colID_tab1, $colID_tab2, \@array_tab1, \@array_tab2);

my @joined_arrays = join_arrays($match_array1_ref,$match_array2_ref);

my $OUT_FHs;

if ($useGzip) {
	# open pipe to Gzip or open text file for writing
	my $out_file = $path2output;
	$out_file =~ s/(.*)\.gz$/$1/;
	my $gz_out_file = $out_file.".gz";
	my $OUT_FHs = new IO::Compress::Gzip ($gz_out_file) or open ">$out_file" or die "Can't open $out_file for writing: $!\n";
} else {
	my $OUT_FHs = open ">$path2output" or die "Can't open $path2output for writing: $!\n";
}

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
	if ( ! -e $path_tab1 ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input file $path_tab1: $!\n"
		);
	}
	if ( ! -e $path_tab2 ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "Cannot find input file $path_tab2: $!\n"
		);
	}
	if (!$path2output ) {
		pod2usage(
			-exitval => 2,
			-verbose => 1,
			-message => "please specify output file name\n"
		);
	}

}

