#!/usr/bin/Rscript

# downlaod genes annotation table from EnsEMBLE BioMart.
# using archived version of a database

# install packages (including BioConductor) if needed
# NOTE: pacman package requires correctly installed CRAN depository address to install required packages.

if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.1.tar.gz",
                                         repos=NULL, method="libcurl")
pacman::p_load( biomaRt )

# use the ensembl mart and the mouse dataset (current version: Ensembl Genes 86)
# ensembl_mm_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
# use archived mart (aug 2014, Ensembl 76)
ensembl_mm_mart76 <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="aug2014.archive.ensembl.org")
#create a filter for all assembled mouse chromosomes
my_chr <- c(1:19, 'M', 'X', 'Y')

# display all attributes
# attributes <- listAttributes(ensembl_mm_mart76)

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
                       mart = ensembl_mm_mart76
)

head(my_annotation)
dim(my_annotation)

# save results as a tab-delimited table
output_dir<-c("/mnt/DB_HDD/Annotation_DB/Mus_musculus.GRCm38/")
filename<-c("mus.musculus.EnsEMBL.GRCm38.annotation.txt")
file<-paste(output_dir, filename, sep="")

write.table(my_annotation.with_ref_seq, file, sep="\t", eol="\n", dec=".", quote=FALSE, row.names = FALSE, col.names = TRUE)


