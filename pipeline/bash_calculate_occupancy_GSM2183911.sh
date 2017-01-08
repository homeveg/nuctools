#!/bin/bash
#
#$ -cwd
#$ -q all.q
#$ -S /bin/bash
#$ -m e 

cd /storage/projects/mESC/GSM2183911_MNase-WT/BED

#Calculate occupancy for each chromosome
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr1.bed --output=1_occup_chr1.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr2.bed --output=1_occup_chr2.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr3.bed --output=1_occup_chr3.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr4.bed --output=1_occup_chr4.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr5.bed --output=1_occup_chr5.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr6.bed --output=1_occup_chr6.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr7.bed --output=1_occup_chr7.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr8.bed --output=1_occup_chr8.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr9.bed --output=1_occup_chr9.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr10.bed --output=1_occup_chr10.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr11.bed --output=1_occup_chr11.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr12.bed --output=1_occup_chr12.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr13.bed --output=1_occup_chr13.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr14.bed --output=1_occup_chr14.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr15.bed --output=1_occup_chr15.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr16.bed --output=1_occup_chr16.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr17.bed --output=1_occup_chr17.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr18.bed --output=1_occup_chr18.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr19.bed --output=1_occup_chr19.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrX.bed --output=1_occup_chrX.occ --window=1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrY.bed --output=1_occup_chrY.occ --window=1