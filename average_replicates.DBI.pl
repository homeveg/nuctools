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

 average_replicates.pl calculates the average occupancy profile and standard deviation based on several replicate occupancy profiles from the working directory and save resulting table, including input occupancy data for individual files. Input *.occ files can be flat or compressed. Resulting extended occupancy file will be saved compressed 

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
use DBI;
use POSIX;

# optional gzip support if modules are installed
my ($ModuleGzipIsLoaded, $ModuleGunzipIsLoaded);
BEGIN { $ModuleGunzipIsLoaded = eval "require IO::Uncompress::Gunzip; 1"; }
BEGIN { $ModuleGzipIsLoaded = eval "require IO::Compress::Gzip; IO::Compress::Gzip->import( qw[gzip] );1"; }

my (%occupancy,@NormFactors);

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
my (@names,@files,@dirs, @all_files,@files_ids);

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
	$file_name =~ s/[\!\"\'\§\$\%\&\/\(\)\=\?\.\;\.\:\-\+\*\#]//g;
	$file_name =~ s/\W//g;
	push(@names, $file_name);
	push(@files_ids, "f$#names");
	}
}
my $file_nr=$#files;


my $dbfile = "$wd/temp.db"; 
my $dsn = "dbi:SQLite:dbname=$dbfile";
my $username     = "";
my $password = "";
my %attr1 = (PrintError=>0, RaiseError=>1, AutoCommit => 0);
my %attr2 = (PrintError=>0, RaiseError=>1, AutoCommit => 1);

# connect without auto-commit
my $dbh = DBI->connect($dsn,$username,$password, \%attr1);

print STDERR "Creating temporary SQLite database...\n";
my $stmt = qq(DROP TABLE IF EXISTS OCCUP );
my $rv = $dbh->do($stmt);
my $sql = join("", "CREATE TABLE OCCUP (COORD INTEGER PRIMARY KEY,", join(" float, ",@files_ids)," float)");
#print STDERR $sql,"\n";
$dbh->do($sql);
$dbh->commit() or die $dbh->errstr;
my $step=100000;

for (my $i=0; $i<=$#files; $i++) {
	my $file_id = $files_ids[$i];
	my $file = $files[$i];
	my $NormFactor = ReadFile("$wd/$file", $file_id, $coordsCol, $occupCol, $dbh, $step, @files_ids);
	push (@NormFactors, $NormFactor);
}

$dbh->commit() or die $dbh->errstr;
$dbh->disconnect;

# re-connect with auto-commit
$dbh = DBI->connect($dsn,$username,$password, \%attr2);

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


my $guery_by_column = $dbh->prepare("SELECT COORD FROM OCCUP");
my $query_by_coord = $dbh->prepare ("SELECT * FROM OCCUP WHERE COORD IS ?");

my $col_ref = $dbh->selectcol_arrayref($guery_by_column);
my @coords = @{ $col_ref};

for (my $i=0; $i<=$#coords; $i++) {
	my $position=$coords[$i];
	$query_by_coord->execute($coords[$i]);
	while (my @row = $query_by_coord->fetchrow_array) {
		( undef, my @temp ) = @row;
		# normalize
		@temp = map { $temp[$_] / $NormFactors[$_] } 0..$#temp;

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

		if ( $i % $step) { }
		else { print STDERR "."; }
	}
}
$dbh->disconnect;

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
    my ($in_file, $filename, $col_coords, $col_occup, $dbh, $step, @names) = @_;
    my $filesize = -s $in_file; #determine file size in bytes
    $filesize = int($filesize/1048576); # filesize in megabytes

    print STDERR "\nReading $in_file file of $filesize MBs.\nPlease wait...";

    #read file with by 4kb chanks
    my $BUFFER_SIZE = 1024*4;
    
	# open compressed occupancy file
	my $inFH;
	if ( $in_file =~ (/.*\.gz$/) ) {
		$inFH = IO::Uncompress::Gunzip->new( $in_file )
		or die "IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip::GunzipError\n";
		
		use constant MIN_COMPRESS_FACTOR => 200;
		my $outer_bytes = -s $in_file;
		my $inner_bytes = get_isize( $in_file );
		$inner_bytes += 4294967296 if( $inner_bytes < $outer_bytes * MIN_COMPRESS_FACTOR );
		$filesize = int($inner_bytes/1048576);
	}
	else { open( $inFH, "<", $in_file ) or die "error: $in_file cannot be opened:$!"; }

    my $buffer = "";
    my $sz_buffer = 0;
    my $timer2 = time();
    # counter for the markers we see
    my $counter = 0;
    
    my $regex_split_newline='\n';
    
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
					
					$sql = "INSERT OR IGNORE INTO OCCUP (COORD) VALUES ($pos);";
					my $add_occup = $dbh->prepare($sql) or die $dbh->errstr;
					$add_occup->execute() or die $add_occup->errstr;
					
					$sql = "UPDATE OCCUP SET $filename=$occup WHERE COORD IS $pos;";			
					$add_occup = $dbh->prepare($sql) or die $dbh->errstr;
					$add_occup->execute() or die $add_occup->errstr;
					
					if ( $counter % $step) { }
					else {  $dbh->commit() or die $dbh->errstr;
						  my $percents = sprintf("%.2f", 100*$offset/(1048576*$filesize));
						  my $processed_fs = sprintf("%.2f", $offset/1048576);
						  my ($progress_seconds,$progress_minutes,$progress_hours,$elapsed_time) = (0,0,0,0);
						  $progress_seconds = time()-$timer2;
						  if ($progress_seconds>=3600) { $progress_hours=floor($progress_seconds/3600); }
						  if ($progress_seconds>=60) { $progress_minutes=floor($progress_seconds/60)-60*$progress_hours; }
						  $progress_seconds=$progress_seconds-3600*$progress_hours-60*$progress_minutes;
						  $elapsed_time=join("", $progress_hours,"h:",$progress_minutes,"m:",$progress_seconds,"s");

						  print STDERR $processed_fs," MBs from $filesize MBs (",$percents,"%) processed in ", $elapsed_time, ".                        \r";}

					push(@all_occups,$occup);
					undef @string;
					$counter++;
        }
        $offset += $n;

        undef @lines;
        $buffer = "";
    }
    
	my ($progress_seconds,$progress_minutes,$progress_hours,$elapsed_time) = (0,0,0,0);
	$progress_seconds = time()-$timer2;
	if ($progress_seconds>=3600) { $progress_hours=floor($progress_seconds/3600); }
	if ($progress_seconds>=60) { $progress_minutes=floor($progress_seconds/60)-60*$progress_hours; }
	$progress_seconds=$progress_seconds-3600*$progress_hours-60*$progress_minutes;
	$elapsed_time=join("", $progress_hours,"h:",$progress_minutes,"m:",$progress_seconds,"s");

    print STDERR "\n",int($offset/1048576), " Mbs processed in ", $elapsed_time, ".\ndone.\n";
	print STDERR "$counter positions added to DB\n";

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

sub get_isize
{
   my ($file) = @_;

   my $isize_len = 4;

   # create a handle we can seek
   my $FH;
   unless( open( $FH, '<:raw', $file ) )
   {
      die "Failed to open $file: $!";
   }
   my $io;
   my $FD = fileno($FH);
   unless( $io = IO::Handle->new_from_fd( $FD, 'r' ) )
   {
      die "Failed to create new IO::Handle for $FD: $!";
   }

   # seek back from EOF
   unless( $io->IO::Seekable::seek( "-$isize_len", 2 ) ) 
   {
      die "Failed to seek $isize_len from EOF: $!"
   }

   # read from here into mod32_isize
   my $mod32_isize;
   unless( my $bytes_read = $io->read( $mod32_isize, $isize_len ) )
   {
      die "Failed to read $isize_len bytes; read $bytes_read bytes instead: $!";
   }

   # convert mod32 to decimal by unpacking value
   my $dec_isize = unpack( 'V', $mod32_isize );

   return $dec_isize;
}