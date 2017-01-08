#!/bin/bash
#$ -cwd
#$ -q all.q
#$ -S /bin/bash
#$ -m e 

mkdir /storage/projects/mESC/100bp_ESC/chr1
mkdir /storage/projects/mESC/100bp_ESC/chr2
mkdir /storage/projects/mESC/100bp_ESC/chr3
mkdir /storage/projects/mESC/100bp_ESC/chr4
mkdir /storage/projects/mESC/100bp_ESC/chr5
mkdir /storage/projects/mESC/100bp_ESC/chr6
mkdir /storage/projects/mESC/100bp_ESC/chr7
mkdir /storage/projects/mESC/100bp_ESC/chr8
mkdir /storage/projects/mESC/100bp_ESC/chr9
mkdir /storage/projects/mESC/100bp_ESC/chr10
mkdir /storage/projects/mESC/100bp_ESC/chr11
mkdir /storage/projects/mESC/100bp_ESC/chr12
mkdir /storage/projects/mESC/100bp_ESC/chr13
mkdir /storage/projects/mESC/100bp_ESC/chr14
mkdir /storage/projects/mESC/100bp_ESC/chr15
mkdir /storage/projects/mESC/100bp_ESC/chr16
mkdir /storage/projects/mESC/100bp_ESC/chr17
mkdir /storage/projects/mESC/100bp_ESC/chr18
mkdir /storage/projects/mESC/100bp_ESC/chr19
mkdir /storage/projects/mESC/100bp_ESC/chr20
mkdir /storage/projects/mESC/100bp_ESC/chr21
mkdir /storage/projects/mESC/100bp_ESC/chr22
mkdir /storage/projects/mESC/100bp_ESC/chrX
mkdir /storage/projects/mESC/100bp_ESC/chrY

cd /storage/projects/mESC/GSM2183911_MNase-WT/BED/
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr1.bed --outdir=/storage/projects/mESC/100bp_ESC/chr1/  --output=chr1_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr2.bed --outdir=/storage/projects/mESC/100bp_ESC/chr2/  --output=chr2_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr3.bed --outdir=/storage/projects/mESC/100bp_ESC/chr3/  --output=chr3_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr4.bed --outdir=/storage/projects/mESC/100bp_ESC/chr4/  --output=chr4_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr5.bed --outdir=/storage/projects/mESC/100bp_ESC/chr5/  --output=chr5_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr6.bed --outdir=/storage/projects/mESC/100bp_ESC/chr6/  --output=chr6_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr7.bed --outdir=/storage/projects/mESC/100bp_ESC/chr7/  --output=chr7_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr8.bed --outdir=/storage/projects/mESC/100bp_ESC/chr8/  --output=chr8_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr9.bed --outdir=/storage/projects/mESC/100bp_ESC/chr9/  --output=chr9_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr10.bed --outdir=/storage/projects/mESC/100bp_ESC/chr10/  --output=chr10_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr11.bed --outdir=/storage/projects/mESC/100bp_ESC/chr11/  --output=chr11_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr12.bed --outdir=/storage/projects/mESC/100bp_ESC/chr12/  --output=chr12_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr13.bed --outdir=/storage/projects/mESC/100bp_ESC/chr13/  --output=chr13_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr14.bed --outdir=/storage/projects/mESC/100bp_ESC/chr14/  --output=chr14_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr15.bed --outdir=/storage/projects/mESC/100bp_ESC/chr15/  --output=chr15_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr16.bed --outdir=/storage/projects/mESC/100bp_ESC/chr16/  --output=chr16_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr17.bed --outdir=/storage/projects/mESC/100bp_ESC/chr17/  --output=chr17_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr18.bed --outdir=/storage/projects/mESC/100bp_ESC/chr18/  --output=chr18_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr19.bed --outdir=/storage/projects/mESC/100bp_ESC/chr19/  --output=chr19_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrX.bed --outdir=/storage/projects/mESC/100bp_ESC/chrX/  --output=chrX_GSM2183911_MNase-WT_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrY.bed --outdir=/storage/projects/mESC/100bp_ESC/chrY/  --output=chrY_GSM2183911_MNase-WT_100bp.occ --window=100

cd /storage/projects/mESC/_West_et_al_2014_nucs/ESC.1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr1.bed --outdir=/storage/projects/mESC/100bp_ESC/chr1/  --output=chr1_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr2.bed --outdir=/storage/projects/mESC/100bp_ESC/chr2/  --output=chr2_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr3.bed --outdir=/storage/projects/mESC/100bp_ESC/chr3/  --output=chr3_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr4.bed --outdir=/storage/projects/mESC/100bp_ESC/chr4/  --output=chr4_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr5.bed --outdir=/storage/projects/mESC/100bp_ESC/chr5/  --output=chr5_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr6.bed --outdir=/storage/projects/mESC/100bp_ESC/chr6/  --output=chr6_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr7.bed --outdir=/storage/projects/mESC/100bp_ESC/chr7/  --output=chr7_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr8.bed --outdir=/storage/projects/mESC/100bp_ESC/chr8/  --output=chr8_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr9.bed --outdir=/storage/projects/mESC/100bp_ESC/chr9/  --output=chr9_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr10.bed --outdir=/storage/projects/mESC/100bp_ESC/chr10/  --output=chr10_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr11.bed --outdir=/storage/projects/mESC/100bp_ESC/chr11/  --output=chr11_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr12.bed --outdir=/storage/projects/mESC/100bp_ESC/chr12/  --output=chr12_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr13.bed --outdir=/storage/projects/mESC/100bp_ESC/chr13/  --output=chr13_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr14.bed --outdir=/storage/projects/mESC/100bp_ESC/chr14/  --output=chr14_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr15.bed --outdir=/storage/projects/mESC/100bp_ESC/chr15/  --output=chr15_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr16.bed --outdir=/storage/projects/mESC/100bp_ESC/chr16/  --output=chr16_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr17.bed --outdir=/storage/projects/mESC/100bp_ESC/chr17/  --output=chr17_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr18.bed --outdir=/storage/projects/mESC/100bp_ESC/chr18/  --output=chr18_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr19.bed --outdir=/storage/projects/mESC/100bp_ESC/chr19/  --output=chr19_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrX.bed --outdir=/storage/projects/mESC/100bp_ESC/chrX/  --output=chrX_West_et_al_2014_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrY.bed --outdir=/storage/projects/mESC/100bp_ESC/chrY/  --output=chrY_West_et_al_2014_rep1_100bp.occ --window=100

cd /storage/projects/mESC/_West_et_al_2014_nucs/ESC.2
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr1.bed --outdir=/storage/projects/mESC/100bp_ESC/chr1/  --output=chr1_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr2.bed --outdir=/storage/projects/mESC/100bp_ESC/chr2/  --output=chr2_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr3.bed --outdir=/storage/projects/mESC/100bp_ESC/chr3/  --output=chr3_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr4.bed --outdir=/storage/projects/mESC/100bp_ESC/chr4/  --output=chr4_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr5.bed --outdir=/storage/projects/mESC/100bp_ESC/chr5/  --output=chr5_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr6.bed --outdir=/storage/projects/mESC/100bp_ESC/chr6/  --output=chr6_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr7.bed --outdir=/storage/projects/mESC/100bp_ESC/chr7/  --output=chr7_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr8.bed --outdir=/storage/projects/mESC/100bp_ESC/chr8/  --output=chr8_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr9.bed --outdir=/storage/projects/mESC/100bp_ESC/chr9/  --output=chr9_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr10.bed --outdir=/storage/projects/mESC/100bp_ESC/chr10/  --output=chr10_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr11.bed --outdir=/storage/projects/mESC/100bp_ESC/chr11/  --output=chr11_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr12.bed --outdir=/storage/projects/mESC/100bp_ESC/chr12/  --output=chr12_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr13.bed --outdir=/storage/projects/mESC/100bp_ESC/chr13/  --output=chr13_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr14.bed --outdir=/storage/projects/mESC/100bp_ESC/chr14/  --output=chr14_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr15.bed --outdir=/storage/projects/mESC/100bp_ESC/chr15/  --output=chr15_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr16.bed --outdir=/storage/projects/mESC/100bp_ESC/chr16/  --output=chr16_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr17.bed --outdir=/storage/projects/mESC/100bp_ESC/chr17/  --output=chr17_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr18.bed --outdir=/storage/projects/mESC/100bp_ESC/chr18/  --output=chr18_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr19.bed --outdir=/storage/projects/mESC/100bp_ESC/chr19/  --output=chr19_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrX.bed --outdir=/storage/projects/mESC/100bp_ESC/chrX/  --output=chrX_West_et_al_2014_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrY.bed --outdir=/storage/projects/mESC/100bp_ESC/chrY/  --output=chrY_West_et_al_2014_rep2_100bp.occ --window=100

cd /storage/projects/mESC/mESC_Zhang2013-nucs/rep1
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr1.bed --outdir=/storage/projects/mESC/100bp_ESC/chr1/  --output=chr1_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr2.bed --outdir=/storage/projects/mESC/100bp_ESC/chr2/  --output=chr2_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr3.bed --outdir=/storage/projects/mESC/100bp_ESC/chr3/  --output=chr3_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr4.bed --outdir=/storage/projects/mESC/100bp_ESC/chr4/  --output=chr4_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr5.bed --outdir=/storage/projects/mESC/100bp_ESC/chr5/  --output=chr5_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr6.bed --outdir=/storage/projects/mESC/100bp_ESC/chr6/  --output=chr6_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr7.bed --outdir=/storage/projects/mESC/100bp_ESC/chr7/  --output=chr7_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr8.bed --outdir=/storage/projects/mESC/100bp_ESC/chr8/  --output=chr8_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr9.bed --outdir=/storage/projects/mESC/100bp_ESC/chr9/  --output=chr9_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr10.bed --outdir=/storage/projects/mESC/100bp_ESC/chr10/  --output=chr10_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr11.bed --outdir=/storage/projects/mESC/100bp_ESC/chr11/  --output=chr11_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr12.bed --outdir=/storage/projects/mESC/100bp_ESC/chr12/  --output=chr12_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr13.bed --outdir=/storage/projects/mESC/100bp_ESC/chr13/  --output=chr13_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr14.bed --outdir=/storage/projects/mESC/100bp_ESC/chr14/  --output=chr14_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr15.bed --outdir=/storage/projects/mESC/100bp_ESC/chr15/  --output=chr15_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr16.bed --outdir=/storage/projects/mESC/100bp_ESC/chr16/  --output=chr16_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr17.bed --outdir=/storage/projects/mESC/100bp_ESC/chr17/  --output=chr17_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr18.bed --outdir=/storage/projects/mESC/100bp_ESC/chr18/  --output=chr18_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr19.bed --outdir=/storage/projects/mESC/100bp_ESC/chr19/  --output=chr19_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrX.bed --outdir=/storage/projects/mESC/100bp_ESC/chrX/  --output=chrX_mESC_Zhang2013_rep1_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrY.bed --outdir=/storage/projects/mESC/100bp_ESC/chrY/  --output=chrY_mESC_Zhang2013_rep1_100bp.occ --window=100


cd /storage/projects/mESC/mESC_Zhang2013-nucs/rep2
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr1.bed --outdir=/storage/projects/mESC/100bp_ESC/chr1/  --output=chr1_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr2.bed --outdir=/storage/projects/mESC/100bp_ESC/chr2/  --output=chr2_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr3.bed --outdir=/storage/projects/mESC/100bp_ESC/chr3/  --output=chr3_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr4.bed --outdir=/storage/projects/mESC/100bp_ESC/chr4/  --output=chr4_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr5.bed --outdir=/storage/projects/mESC/100bp_ESC/chr5/  --output=chr5_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr6.bed --outdir=/storage/projects/mESC/100bp_ESC/chr6/  --output=chr6_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr7.bed --outdir=/storage/projects/mESC/100bp_ESC/chr7/  --output=chr7_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr8.bed --outdir=/storage/projects/mESC/100bp_ESC/chr8/  --output=chr8_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr9.bed --outdir=/storage/projects/mESC/100bp_ESC/chr9/  --output=chr9_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr10.bed --outdir=/storage/projects/mESC/100bp_ESC/chr10/  --output=chr10_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr11.bed --outdir=/storage/projects/mESC/100bp_ESC/chr11/  --output=chr11_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr12.bed --outdir=/storage/projects/mESC/100bp_ESC/chr12/  --output=chr12_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr13.bed --outdir=/storage/projects/mESC/100bp_ESC/chr13/  --output=chr13_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr14.bed --outdir=/storage/projects/mESC/100bp_ESC/chr14/  --output=chr14_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr15.bed --outdir=/storage/projects/mESC/100bp_ESC/chr15/  --output=chr15_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr16.bed --outdir=/storage/projects/mESC/100bp_ESC/chr16/  --output=chr16_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr17.bed --outdir=/storage/projects/mESC/100bp_ESC/chr17/  --output=chr17_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr18.bed --outdir=/storage/projects/mESC/100bp_ESC/chr18/  --output=chr18_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chr19.bed --outdir=/storage/projects/mESC/100bp_ESC/chr19/  --output=chr19_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrX.bed --outdir=/storage/projects/mESC/100bp_ESC/chrX/  --output=chrX_mESC_Zhang2013_rep2_100bp.occ --window=100
perl /storage/projects/nuctools.2.0/bed2occupancy_average.pl --input=chrY.bed --outdir=/storage/projects/mESC/100bp_ESC/chrY/  --output=chrY_mESC_Zhang2013_rep2_100bp.occ --window=100