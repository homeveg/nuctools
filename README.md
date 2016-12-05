# NucTools

## NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data
 
Recent advancements in high-throughput sequencing methods create a vast amount of data in which numerous chromatin features are mapped along the genome. The results are frequently analysed by creating binarized data sets that link the presence/absence of a given readout to specific genomic loci. It is currently a challenge in the field to cope with continuous distributions of deep sequencing chromatin readouts, and to integrate the different types of datasets to reveal linkages between them. Here we introduce the NucTools suite of Perl scripts and a stand-alone visualisation program for a nucleosome-centred downstream analysis of deep sequencing data. NucTools accounts for the continuous distribution of nucleosome occupancy. Furthermore, it is useful to associate nucleosome occupancy with other chromatin features like transcription factor binding, histone modifications or DNA methylation.

-------------------------------------------------------------------
## SYSTEM REQUIREMENT

Linux (2.6 kernel or later) or Mac (OSX 10.6 Snow Leopard or later) operating system
with minimal 16 GB of RAM is recommended*. Perl v5.8 or above is required.

The C/C++ compiling enviroment might be required for installing dependencies, such as bedtools. Systems may vary.
Please assure that your system has the essential software building packages (e.g. build-essential
for Fedora, XCODE for Mac...etc) installed properly before running the installing 
script.

NucTools was tested successfully on our Linux servers (CentOS release 6.7 w/ Perl v5.10.1; 
Fedora release 22 w/ Perl v5.20.3) and Macbook Pro laptops (MAC OSX 10.11 w/ XCODE v5.1).

*Memory requirments depends on the experimental system. For big genomes performance will increase greatly on machines with more memory. For example, processing of 22 mouse chromosomes, with the sequencing library size about 75 000 000 reads occupy at peak load around 60-70 Gb of RAM. Important to mention, the performance is very dependent on HDD read-write speed. Therefore the parallel running of many samples at once recomended only on the server-like system or computational clusters with RAID arrays, allowing real multithreaded read-write access to HDDs array.

-------------------------------------------------------------------
## QUICK START 

This is an example of profiling a "test.bed" file using NucTools. The testing BED file comes along with the NucTools package
in the "test" directory. More details are stated in the INSTRUCTION section.

1. Obtaining NucTools package:
   
        $ git clone https://github.com/homeveg/nuctools.git NucTools
   
2. Installing NucTools:
        the NucTools package does not require installation. It is a collection of individual Perl scripts which can be executed individually.

3. Generate a genome annotation table using provided R script:
   
        $ Rscript misc/LoadAnnotation.BioMart.R

4. Prepare BED files from BAM files with external application (optionally) 
    1. merge multiple replicates to one BAM file and sort by reads name:
   
        $ samtools merge -n /Path_to_folder_with/BAM/test_sorted.bam /Path_to_folder_with/BAM/test.rp1.bam \ 
        /Path_to_folder_with/BAM/test.rp2.bam /Path_to_folder_with/BAM/test.rp3.bam
        
    2. converted sorted BAM file to BED using bowtie2bed.pl script or with external app bedtools:
   
        $ perl -w bowtie2bed.pl -i /Path_to_folder_with/BAM/test_sorted.bam -verbose > /Path_to_folder_with/BED/test_sorted.bed.gz
        $ bedtools bamtobed -i /Path_to_folder_with/BAM/test_sorted.bam | pigz > /Path_to_folder_with/BED/test_sorted.bed.gz

5. Running NucTools:
    1. Extend single-end reads to the average DNA fragment size
   
        $ extend_SE_reads.pl -in test.bed -out test.ext.bed.gz -fL 147
        
    2. Extract individual chromosomes from whole genome BED
   
        $ extract_chr_bed.pl -in test.ext.bed.gz -out test -d /Path_to_folder_with/BED/ -p chr 
        
    3. Convert all BED files to occupancy OCC files averaging nucleosomes occupancy values over the window of size 10
   
        $ bed2occupancy_average.pl -in /Path_to_folder_with/BED/ -odir /Path_to_folder_with/OCC -dir -use -w 10
        
    4. Calculate aggregate profiles and aligned occupancy matrices for each chromosome individually
    
        $ aggregate_profile.pl -reg genome_annotation.tab -idC 0 -chrC 4 --strC 7 -sC 8 -eC 9 -pbN -lsN -lS <SeqLibSize> \
        -chr 1 -al /Path_to_folder_with/OCC/chr1.test.occ_matrix -av /Path_to_folder_with/OCC/chr1.test.aggregate \
        -in /Path_to_folder_with/OCC/chr1.test.w10.occ.gz -upD 1000 -downD 1000
        
    5. Paste together aggregate profiles of each chromosome in one file and add a header
       * Generate a header
   
        $ ls /Path_to_folder_with/OCC/*1000_1000.txt | perl -n -e 'if(/.*(chr.*)\.test.*/gm) { print $1, "\t"; }' | \
        perl -n -e 'if( /(.*)\t$/g )  { print $1}' > /Path_to_folder_with/OCC/test.all.occ.txt
        $ echo "" >> /Path_to_folder_with/OCC/test.all.occ.txt
        
       * paste all aggregate profiles to one file
   
        $ paste /Path_to_folder_with/OCC/*1000_1000.txt >> /Path_to_folder_with/OCC/test.all.occ.txt

6. Optionally: 
Visualize aggregate profiles and run K-mean cluster analysis on aligned occupancy matrixes with the MatLab-based ClusterMaps Building Tool (provided as a part of NucTools package). 
[Download link](https://github.com/homeveg/nuctools/raw/master/CMB/ClusterMaps_Builder_src.zip)


-------------------------------------------------------------------
### Installation

the NucTools suite for a nucleosome-centred downstream analysis of deep sequencing data is primarily Perl-based, and require at least Perl v5.8 with dependencies installed properly (listed in README_FULL.md).
A visualisation program is written on MatLab and requires either full MatLab insatllation or can be provided as a standolone application with web-installer compiled for MacOS X 10.9 or later.
NucTools utilize whole genome [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) files.

**Optional external applications:**
    
- [SamTools](http://samtools.sourceforge.net/) - merge, sort and convert BAM files
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html) - convert BAM to BED
- [PIGZ](http://zlib.net/pigz/) - a parallel implementation of gzip for modern multi-processor, multi-core machines

-------------------------------------------------------------------
### Running NucTools

One can divide the NucTools pipeline to 3 major steps: (1) prepare input OCC files (2) calculate aggregate profiles and alligned occupancy matrixes (3) followup analysis and results visualization. 
For the moment we don't have a wrapper to run all 3 steps automatically so, each step should be executed separately and, in turn, consists of several intermediate steps.

We will provide a testing BAM file, command line bash file, and example output in ./test. Below is a set of commands that runs "test.bam" through NucTools pipeline:

        $ samtools sort -n ./test/test_sorted.bam ./test/test.bam
        $ bowtie2bed.pl -i ./test/test_sorted.bam --verbose > ./test/test_sorted.bed.gz
        $ extend_SE_reads.pl -in ./test/test.bed -out ./test/test.ext.bed.gz
        $ extract_chr_bed.pl -in ./test/test.ext.bed.gz -out test/BED -d ./test -p chr 
        $ bed2occupancy_average.pl -in ./test/BED -odir ./test/OCC -dir -use -w 10
        $ aggregate_profile.pl -reg genome_annotation.txt -idC 0 -chrC 4 -strC 7 -sC 8 -eC 9 -pbN -lsN -lS 75000000 -chr 1 -al ./test/OCC/chr1.test.occ_matrix -av ./test/OCC/chr1.test.aggregate -in ./test/OCC/chr1.test.w10.occ.gz -upD 1000 -downD 1000

In the example above test.bam file is sorted by the reads names and converted to test_sorted.bed.gz file. As long as in the example case we are dealing with single-end ilumina sequencing reads, with expected read length of 100, we extending each read to the length of 100 using extend_SE_reads.pl. Resulting whole genome BED file is splitted to chromosomes with extract_chr_bed.pl script and all per-chromosome bed files are converted to OCC files with bed2occupancy_average.pl, using running window 10 (use -w 0 for bep-base resolution). Last step is to generate an aggregate profiles using aggregate_profile.pl script 

-------------------------------------------------------------------
### Interpreting Results

NucTools bed2occupancy_average.pl scritp generates 2 types of the output. The file containing one column of numbers, corresponding to an aggregate profile (could be visualized, for example in Excel) and the tab-delimited text file containing all occupancy data for each region of interest (or transcript), aligned at annotated strat. These matrixes can be later visualized with our ClusterMaps Building Tool.


-------------------------------------------------------------------
## NucTools scripts

### Initial data transformation

   * ### bowtie2bed.pl
   takes as an input standard SAM, BAM or MAP file and converts to the gzip-compressed BED file. The programm require samtools installed in PATH to be able to work with BAM files
   
        $ perl -w bowtie2bed.pl --input=accepte_hits.bam --output=sample.bed.gz [--verbose --help]
   
   * ### extend_SE_reads.pl
   Extends single-end reads by the user-defined value of the average DNA fragment length. Script works with compressed or uncompressed BED files and save output as compress *.BED.GZ
   
        $ perl -w extend_SE_reads.pl -in <in.bed> -out <out.bed> -fL <fragment length> \
        [-cC <column Nr.> -sC <column Nr.> -eC <column Nr.> -strC <column Nr.> ] [--help] 
   
   * ### extend_PE_reads.pl
   Takes as an input BED file with mapped paired-end reads (two lines per paired read) sorted according to the read name and reformat it by creating a smaller BED file with one line per nucleosome in the following format: (1) chromosome, (2) nucleosome start, (3) nucleosome end, (4) nucleosome length
    
        $ perl -w extend_PE_reads.pl -in <in.bed> -out <out.bed> [--help] 
   
   * ### calc_fragment_length.pl
   Estimates mean fragment length for a single-emd sequencing (can be used for ssingle end reads extention)
    
        $ perl -w perl -w calc_fragment_length.pl --input=<in.bed> --output=<filtered.txt> [--delta=<N> --apply_filter \
        --filtering_threshold=<N> --pile=<N> --fix_pile_size ] [--chromosome_col=<column Nr.> --start_col=<column Nr.> \
        --end_col=<column Nr.> --strand_col=<column Nr.> --help]
   
   * ### extract_chr_bed.pl
   Splits whole genome BED file with mapped reads into smaller BED files per each chromosome
    
        $ perl -w extract_chr_bed.pl -in all_data.bed.gz -out output_name_template -p [<pattern>] [--help] 
   
   * ### bed2occupancy_average.pl
   Calculates genome-wide nucleosome occupancy, based on the BED file with sequencing reads. It converts BED files for all or specified chromosomes. The running window occupancy file (*.OCC) is a text file containing normalized reads frequency distribution along each chromosome in the running window.
    
        $ perl -w bed2occupancy_average.pl --input=<in.bed.gz> --output=<out.occ.gz> \
        [--outdir=<DIR_WITH_OCC> --chromosome_col=<column Nr.> --start_col=<column Nr.> --end_col=<column Nr.> \
        --strand_col=<column Nr.> --window=<running window size> --consider_strand --ConvertAllInDir --help]
   

### Core scripts

   * ### aggregate_profile.pl
   Calculates aggregate profile of sequencing read density around genomic regions. As an input it utilze a tab-delimited text file or BED file with coordinates of genomic features (promoters, enhancers, chromatin domains, TF binding sites, etc), and the OCC files with continuous chromosome-wide occupancy (nucleosome occupancy, TF distribution, etc). Calculates normalized occupancy profiles for each of the features, as well as the aggregate profile representing the average occupancy centred at the middle of the feature
    
        $ perl -w aggregate_profile.pl --input=<in.occ.gz> --regions=<annotations.txt> [--expression=<gene_expression.rpkm>] \ 
        --aligned=<output.aligned.tab.gz> --average_aligned=<output.aggregare.txt> \ 
        [--path2log=<AggregateProfile.log> --region_start_column=<column Nr.> --region_end_column=<column Nr.> \
        --strand_column=<column Nr.> --chromosome_col=<column Nr.> --GeneId_column=<column Nr.> \
        --Expression_columnID=<column Nr.> --Methylation_columnID=<column Nr.> --Methylation_columnID2=<column Nr.> \
        --upstream_delta=<column Nr.> --downstream_delta==<column Nr.> --upper_threshold=<column Nr.> --lower_threshold=<column Nr.> \
        --Methylation_threshold=<value|range_start-range_end> --overlap=<length> --library_size=<Nr.> \
        --remove_minus_strand | --ignore_strand | --fixed_strand=[plus|minus] --invert_strand --input_occ --score --dont_save_aligned \
        --Cut_tail --chromosome=chrN --AgregateProfile --GeneLengthNorm --LibsizeNorm --PerBaseNorm --useCentre \
        --use_default --verbose --help ]

   * ### average_replicates.pl
   Calculates the average occupancy profile and standard deviation based on several replicate occupancy profiles from the working directory and save resulting table, including input occupancy data for individual files. Input *.occ files can be flat or compressed. Resulting extended occupancy file will be saved compressed
    
        $ perl -w average_replicates.pl --dir=<path to working dir> --output=<path to results file> --coordsCol=0 \
        --occupCol=1 --pattern="occ.gz" --printData --sum [--help]  
   
   * ### calc_fragment_length.pl
   Estimates mean fragment length for a single-emd sequencing library
    
        $ perl -w calc_fragment_length.pl --input=<in.bed> --output=<filtered.txt> \
        [--delta=<N> --apply_filter --filtering_threshold=<N> --pile=<N> --fix_pile_size ] \ 
        [--chromosome_col=<column Nr.> --start_col=<column Nr.> --end_col=<column Nr.> --strand_col=<column Nr.> --help]
   
   * ### nucleosome_repeat_length.pl
   Calculates frequency of nucleosome-nucleosome distances to determine the nucleosome repeat length
     
        $ perl -w nucleosome_repeat_length.pl --input=<in.bed> --output=<filtered.txt> \
        [--delta=<N> --apply_filter --filtering_threshold=<N> --pile=<N> --fix_pile_size ] \
        [--chromosome_col=<column Nr.> --start_col=<column Nr.> --end_col=<column Nr.> --strand_col=<column Nr.> --help]
   
   * ### stable_nucs_replicates.pl
   Finds stable and fussy nucleosomes using all replicates for the same experimental condition
     
        $ perl -w stable_nucs_replicates.pl --input=<path to input DIR> --output=<out.bed> --chromosome=chr1 \
        [-coordsCol=0 -occupCol=2 -StableThreshold=0.5 --printData ] [--help] 
   
   
### Vizualization and additional scripts

   * ### LoadAnnotation.BioMart.R
   R script to retrieve genes annotation from EnsEMBL using Bioconductor BioMart package. Genes annotation table, particulary TSS/TTS coordinates, chromosomes and strand inforamtion is used with aggregate_profile.pl as a genomic features table.
   
   * ### NRL.loess.estimation.R
   Peak detection R script to estimate NRL based on nucleosome_repeat_length.pl output.
   
   * ### CMB - Cluster Maps Builder
   Aggregate profile and aligned occupancy matrix visualizer. MatLab-based stand-alown GUI application, compiled to run on MacOS X

-------------------------------------------------------------------
## Additional information

Additional information, publications references and short description of each script from the toolbox can be found here:
- http://www.generegulation.info/index.php/nuctools (external link)
- https://homeveg.github.io/nuctools/ (GitHub pages)

**PLEASE NOTE**:
*we are currently preparing a manuscript and will add here more information, usage instruction and test data set soon.
The project is under development now.*

### Planned modifications

There are several changes has been planned and partially implemented:

- implement GZIP support for all input/output files generated by the pipeline scripts
- implement gapped OCC format support to all relevant scripts
- extend NucTools documentation and usage examples with Perl POD package
- unify pipeline scripts input parameters using Perl GETOPT package

### Future possible modifications

- BAM files support
- parallel processing (beautiful codes snippets for implementation of parallel processing of BAM files with Perl one can find here: https://genomebytes.wordpress.com/2013/07/24/multi-thread-access-of-bam-files-using-perl-and-samtools/ )
- NucTools automated package installation with make

## Major changes (in porogres):
- GZIP support implementation - all input files can be compressed and all newly created output files saved as *.gz
- Handle gapped occupancy files: to save space and memory remove all regions/bases with 0 counts
- implement POD documentation (with Pod::Usage)
- command line parametrs allow long and short options (with Getopt::Long::GetOptions)

### Developers: 
Yevhen Vainshtein and Vladimir B. Teif
