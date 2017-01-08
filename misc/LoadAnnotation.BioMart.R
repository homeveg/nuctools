#!/usr/bin/env Rscript
# download genes annotation table from EnsEMBLE BioMart.
# using archivedor current version of a database

# load required packages
if (!require("devtools")) install.packages("https://cran.r-project.org/src/contrib/devtools_1.12.0.tar.gz",
                                           repos=NULL, method="libcurl")
if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.1.tar.gz",
                                         repos=NULL, method="libcurl")
pacman::p_load( optparse, biomaRt )

# read command line arguments
option_list = list(
  make_option("--biomart", type="character", default="ENSEMBL_MART_ENSEMBL", 
              help="select biomaRT database to use [default= %default]", metavar="bimaRT type"),
  make_option("--dataset", type="character", default="mmusculus_gene_ensembl", 
              help="biomaRT dataset. Use listDatasets(ensembl) command to see full list of all Ensembl species [default= %default]", 
              metavar="bimaRT dataset"),
  make_option("--host", type="character", default="aug2014.archive.ensembl.org", 
              help="Accessing archives through specifying the archive host [default= %default]. To retrieve data fromcurrent EnsEmbl release use host=current", 
              metavar="archive host addres"),
  make_option("--annotation", type="character", default="genes_annotation.txt", 
              help="resulting genes annotation table [default= %default]", metavar="file name"),
  make_option("--dir", type="character", default=NULL, 
              help="path to output directory", metavar="path"),
  make_option("--ChrNrs", type="integer", default=19, 
              help="Number of chromosomes to retrieve except X,Y and M [default= %default]", metavar="chromosomes")
); 
""
# parse arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check if required arguments specified
if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Please specify required parameter: path to a output directory.", call.=FALSE)
}

if (opt$host == "current") {
  ensembl_mart <- useMart(opt$biomart, dataset=opt$dataset)
} else {
  ensembl_mart <- useMart(opt$biomart, dataset=opt$dataset, host=opt$host)
}

#create a filter for all assembled chromosomes
my_chr <- c(1:opt$ChrNrs, 'M', 'X', 'Y')

# display all attributes
# attributes <- listAttributes(ensembl_mart)

# list RefSeq-related attributes
# attributes[grep(pattern="refseq", x=attributes$description, ignore.case=T),]

# list RefSeq-related attributes
# attributes[grep(pattern="ensembl", x=attributes$description, ignore.case=T),]

# requested attributes
my_attributes <- c('refseq_mrna', 'external_gene_name', 'ensembl_gene_id', 'ensembl_transcript_id', 
                   'chromosome_name', 'start_position', 'end_position', 'strand', 
                   'transcript_start', 'transcript_end')

# request all attributes listed in my_attributes from archived Ensembl 76 database and filter for all assembled mouse chromosomes
my_annotation <- getBM(attributes = my_attributes,
                       filters = 'chromosome_name',
                       values = my_chr,
                       mart = ensembl_mart
)

head(my_annotation)
dim(my_annotation)

# save results as a tab-delimited table
output_dir<-opt$dir
filename<-opt$annotation
file<-paste(output_dir, filename, sep="")

write.table(my_annotation, file, sep="\t", eol="\n", dec=".", quote=FALSE, row.names = FALSE, col.names = TRUE)

print(paste("Annotation file saved to",file))


