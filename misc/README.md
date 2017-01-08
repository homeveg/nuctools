# NucTools R scripts

## Disclaimer

The following collection of R scripts is part of NucTools 2.0 release. At the moment it is at constant development and therefore may have some instability and bugs. If you find new bugs please contact yevhen.vainshtein@igb.fraunhofer.de

-------------------------------------------------------------------

## SYSTEM REQUIREMENT

- R statistical language (https://www.r-project.org/)
Required packages:
Each script provides an automated packages-installation pipeline based on R package “pacman” and normally does not require manual installation of R packages. 
In some systems one might need to install the “libcurl” and/or “devtools” packages manually. Installation instructuoin can be found here:
-	libcurl: https://cran.r-project.org/web/packages/curl/index.html
-	devtools: https://cran.r-project.org/web/packages/devtools/README.html

-------------------------------------------------------------------

## Introduction

The collection consists of the following R scripts:
- LoadAnnotation.BioMart.R
- match_2tables_byID.R
- plotNRL.R

All scripts can be executed from the command line with parameters specifying input/output files and other run options. Script LoadAnnotation.BioMart.R may require some basic knowledge of R programming language to introduce changes to the script code or list possibly databases and species.

-------------------------------------------------------------------

## Quick start

From the terminal window call script using following command:

    $ Rscript <script name.R> --help

### LoadAnnotation.BioMart.R

The script LoadAnnotation.BioMart.R allows downloading a detailed genes annotation table from EnsEMBL Biomart database (www.ensembl.org/biomart)
The script is using Bioconductor package boimaRT to retrieve genes annotation from latest or archived databases (https://bioconductor.org/packages/release/bioc/html/biomaRt.html) 
By default script LoadAnnotation.BioMart.R will download mouse genes annotation from archived Ensembl 76 database and save it to a user specified directory.

To start working with the script use the following command:

    $ Rscript LoadAnnotation.BioMart.R --help

### match_2tables_byID.R

ClusterMapsBuilder (CMB) has an option to apply the sorting order after K-mean clustering from one dataset to another. The prerequisite of such analysis is two identically sorted aligned occupancy matrices of the same length.
The script match_2tables_byID.R is developed to prepare aligned occupancy matrices (the result of aggregate_profile.pl script) for further analysis with CMB. It resorting each of 2 input tables by the first column, containing unique ID and saving to the output only lines with matching IDs.

Execute the following command to get detailed usage instruction:

    $ Rscript match_2tables_byID.R --help

### plotNRL.R

The distribution of nucleosome start-to-start distances determined by nucleosome_repeat_length.pl can be analysed by an R script plotNRL.R, which extracts peak coordinates and performs linear fitting; the slope of the line gives the NRL.

Execute the following command to get detailed usage instruction:

    $ Rscript plotNRL.R --help

### Developers: 
Yevhen Vainshtein and Vladimir B. Teif
