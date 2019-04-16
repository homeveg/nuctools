#!/usr/bin/env Rscript

# load required packages
if (!require("devtools")) install.packages("https://cran.r-project.org/src/contrib/devtools_2.0.2.tar.gz",
                                           repos=NULL, method="libcurl")
if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.5.1.tar.gz",
                                         repos=NULL, method="libcurl")
pacman::p_load( optparse )

##############################################################################

############################
# read command line options
############################

# read command line arguments
option_list = list(
  make_option("--input1", type="file", default=NULL, 
              help="Occupancy matrix 1", metavar="file name"),
  make_option("--input2", type="character", default=NULL, 
              help="Occupancy matrix 2", metavar="file name"),
  make_option("--output1", type="character", default="out1.occ.txt", 
              help="Sorted and matched occupancy matrix 1", metavar="file name"),
  make_option("--output2", type="character", default="out2.occ.txt", 
              help="Sorted and matched occupancy matrix 2", metavar="file name"),
  make_option("--dir", type="character", default=NULL, 
              help="path to working directory", metavar="path")
); 

# parse arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check if required arguments specified
if (is.null(opt$input1)){
  print_help(opt_parser)
  stop("Please specify input occupancy matrix 1 file name.", call.=FALSE)
}

if (is.null(opt$input2)){
  print_help(opt_parser)
  stop("Please specify input occupancy matrix 2 file name.", call.=FALSE)
}

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Please specify path to a working directory.", call.=FALSE)
}

wd <- opt$dir


#TSS
in_files <- c(paste(wd,opt$input1,sep="/"),paste(wd,opt$input2,sep="/") )
out_files <- c(  paste(wd,opt$output1,sep="/"),paste(wd,opt$output2,sep="/") )

##########################
# compare translatome pipelines
##########################

data1<-as.matrix(read.table(in_files[1], header=FALSE, sep="\t", dec=".", fill=TRUE))
data2<-as.matrix(read.table(in_files[2], header=FALSE, sep="\t", dec=".", fill=TRUE))

rownames(data1)<-data1[,1]
rownames(data2)<-data2[,1]

rownames_data1<-rownames(data1)
rownames_data2<-rownames(data2)

intersected_names<-intersect(rownames_data1,rownames_data2)

common<-data1[match(intersected_names,rownames_data1),]
common_notNA<-common[!is.na(common[,3]),]
common_notNA<-common_notNA[,-1]
common_notNA[is.na(common_notNA)]<-0
data1_new<-common_notNA

common<-data2[match(intersected_names,rownames_data2),]
common_notNA<-common[!is.na(common[,3]),]
common_notNA<-common_notNA[,-1]
common_notNA[is.na(common_notNA)]<-0
data2_new<-common_notNA

write.table(data1_new, file=out_files[1], quote=FALSE, eol="\n", na = "0",dec=".", sep="\t", row.names=TRUE, col.names =FALSE)
write.table(data2_new, file=out_files[2], quote=FALSE, eol="\n", na = "0", dec=".", sep="\t", row.names=TRUE, col.names =FALSE)
