# NucTools

## NucTools is a software package for analysis of chromatin feature occupancy profiles from high-throughput sequencing data
 
Recent advancements in high-throughput sequencing methods create a vast amount of data in which numerous chromatin features are mapped along the genome. The results are frequently analysed by creating binarized data sets that link the presence/absence of a given readout to specific genomic loci. It is currently a challenge in the field to cope with continuous distributions of deep sequencing chromatin readouts, and to integrate the different types of datasets to reveal linkages between them. Here we introduce the NucTools suite of Perl scripts and a stand-alone visualisation program for a nucleosome-centred downstream analysis of deep sequencing data. NucTools accounts for the continuous distribution of nucleosome occupancy. Furthermore, it is useful to associate nucleosome occupancy with other chromatin features like transcription factor binding, histone modifications or DNA methylation.

## Additional information

Additional information and short description of each script from the toolbox can be found here:
http://www.generegulation.info/index.php/nuctools

## PLEASE NOTE:
we are currently preparing manuscript and will add here more information, usage instruction and test data set soon.
The project is under development now.

## Major changes (in porogres):
- GZIP support implementation - all input files can be compressed and all newly created output files saved as *.gz
- Handle gapped occupancy files: to save space and memory remove all regions/bases with 0 counts
- implement POD documentation (with Pod::Usage)
- command line parametrs allow long and short options (with Getopt::Long::GetOptions)

### Developers: 
Yevhen Vainshtein and Vladimir B. Teif
