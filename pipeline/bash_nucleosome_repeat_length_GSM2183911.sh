#!/bin/bash
#
#$ -cwd
#$ -q all.q
#$ -S /bin/bash
#$ -m e 

cd /storage/projects/mESC/GSM2183911_MNase-WT/BED

perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr1.bed --output=nuc-nuc_v5_chr1.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr2.bed --output=nuc-nuc_v5_chr2.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr3.bed --output=nuc-nuc_v5_chr3.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr4.bed --output=nuc-nuc_v5_chr4.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr5.bed --output=nuc-nuc_v5_chr5.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr6.bed --output=nuc-nuc_v5_chr6.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr7.bed --output=nuc-nuc_v5_chr7.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr8.bed --output=nuc-nuc_v5_chr8.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr9.bed --output=nuc-nuc_v5_chr9.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr10.bed --output=nuc-nuc_v5_chr10.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr11.bed --output=nuc-nuc_v5_chr11.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr12.bed --output=nuc-nuc_v5_chr12.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr13.bed --output=nuc-nuc_v5_chr13.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr14.bed --output=nuc-nuc_v5_chr14.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr15.bed --output=nuc-nuc_v5_chr15.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr16.bed --output=nuc-nuc_v5_chr16.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr17.bed --output=nuc-nuc_v5_chr17.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr18.bed --output=nuc-nuc_v5_chr18.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chr19.bed --output=nuc-nuc_v5_chr19.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chrX.bed --output=nuc-nuc_v5_chrX.txt --delta=2000 --MaxNr=10000000 --apply_filter
perl /storage/projects/nuctools.2.0/nucleosome_repeat_length.pl --input=chrY.bed --output=nuc-nuc_v5_chrY.txt --delta=2000 --MaxNr=10000000 --apply_filter