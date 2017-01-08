#!/bin/bash
#$ -cwd
#$ -q all.q
#$ -S /bin/bash
#$ -m e

cd /storage/projects/mESC/100bp_MEF_vs_ESC

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr1 --output=chr1_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr1"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr1 --output=chr1_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr1"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr1_average_ESC_100bp.txt --input2=chr1_average_MEF_100bp.txt --output1=chr1_MEF_less_ESC_0.99.txt --output2=chr1_MEF_more_ESC_0.99.txt --chromosome=chr1 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr2 --output=chr2_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr2"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr2 --output=chr2_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr2"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr2_average_ESC_100bp.txt --input2=chr2_average_MEF_100bp.txt --output1=chr2_MEF_less_ESC_0.99.txt --output2=chr2_MEF_more_ESC_0.99.txt --chromosome=chr2 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr3 --output=chr3_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr3"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr3 --output=chr3_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr3"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr3_average_ESC_100bp.txt --input2=chr3_average_MEF_100bp.txt --output1=chr3_MEF_less_ESC_0.99.txt --output2=chr3_MEF_more_ESC_0.99.txt --chromosome=chr3 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr4 --output=chr4_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr4"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr4 --output=chr4_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr4"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr4_average_ESC_100bp.txt --input2=chr4_average_MEF_100bp.txt --output1=chr4_MEF_less_ESC_0.99.txt --output2=chr4_MEF_more_ESC_0.99.txt --chromosome=chr4 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr5 --output=chr5_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr5"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr5 --output=chr5_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr5"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr5_average_ESC_100bp.txt --input2=chr5_average_MEF_100bp.txt --output1=chr5_MEF_less_ESC_0.99.txt --output2=chr5_MEF_more_ESC_0.99.txt --chromosome=chr5 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr6 --output=chr6_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr6"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr6 --output=chr6_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr6"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr6_average_ESC_100bp.txt --input2=chr6_average_MEF_100bp.txt --output1=chr6_MEF_less_ESC_0.99.txt --output2=chr6_MEF_more_ESC_0.99.txt --chromosome=chr6 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr7 --output=chr7_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr7"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr7 --output=chr7_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr7"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr7_average_ESC_100bp.txt --input2=chr7_average_MEF_100bp.txt --output1=chr7_MEF_less_ESC_0.99.txt --output2=chr7_MEF_more_ESC_0.99.txt --chromosome=chr7 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr8 --output=chr8_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr8"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr8 --output=chr8_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr8"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr8_average_ESC_100bp.txt --input2=chr8_average_MEF_100bp.txt --output1=chr8_MEF_less_ESC_0.99.txt --output2=chr8_MEF_more_ESC_0.99.txt --chromosome=chr8 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr9 --output=chr9_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr9"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr9 --output=chr9_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr9"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr9_average_ESC_100bp.txt --input2=chr9_average_MEF_100bp.txt --output1=chr9_MEF_less_ESC_0.99.txt --output2=chr9_MEF_more_ESC_0.99.txt --chromosome=chr9 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr10 --output=chr10_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr10"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr10 --output=chr10_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr10"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr10_average_ESC_100bp.txt --input2=chr10_average_MEF_100bp.txt --output1=chr10_MEF_less_ESC_0.99.txt --output2=chr10_MEF_more_ESC_0.99.txt --chromosome=chr10 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr11 --output=chr11_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr11"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr11 --output=chr11_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr11"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr11_average_ESC_100bp.txt --input2=chr11_average_MEF_100bp.txt --output1=chr11_MEF_less_ESC_0.99.txt --output2=chr11_MEF_more_ESC_0.99.txt --chromosome=chr11 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr12 --output=chr12_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr12"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr12 --output=chr12_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr12"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr12_average_ESC_100bp.txt --input2=chr12_average_MEF_100bp.txt --output1=chr12_MEF_less_ESC_0.99.txt --output2=chr12_MEF_more_ESC_0.99.txt --chromosome=chr12 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr13 --output=chr13_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr13"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr13 --output=chr13_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr13"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr13_average_ESC_100bp.txt --input2=chr13_average_MEF_100bp.txt --output1=chr13_MEF_less_ESC_0.99.txt --output2=chr13_MEF_more_ESC_0.99.txt --chromosome=chr13 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr14 --output=chr14_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr14"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr14 --output=chr14_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr14"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr14_average_ESC_100bp.txt --input2=chr14_average_MEF_100bp.txt --output1=chr14_MEF_less_ESC_0.99.txt --output2=chr14_MEF_more_ESC_0.99.txt --chromosome=chr14 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr15 --output=chr15_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr15"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr15 --output=chr15_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr15"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr15_average_ESC_100bp.txt --input2=chr15_average_MEF_100bp.txt --output1=chr15_MEF_less_ESC_0.99.txt --output2=chr15_MEF_more_ESC_0.99.txt --chromosome=chr15 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr16 --output=chr16_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr16"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr16 --output=chr16_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr16"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr16_average_ESC_100bp.txt --input2=chr16_average_MEF_100bp.txt --output1=chr16_MEF_less_ESC_0.99.txt --output2=chr16_MEF_more_ESC_0.99.txt --chromosome=chr16 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr17 --output=chr17_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr17"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr17 --output=chr17_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr17"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr17_average_ESC_100bp.txt --input2=chr17_average_MEF_100bp.txt --output1=chr17_MEF_less_ESC_0.99.txt --output2=chr17_MEF_more_ESC_0.99.txt --chromosome=chr17 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr18 --output=chr18_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr18"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr18 --output=chr18_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr18"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr18_average_ESC_100bp.txt --input2=chr18_average_MEF_100bp.txt --output1=chr18_MEF_less_ESC_0.99.txt --output2=chr18_MEF_more_ESC_0.99.txt --chromosome=chr18 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chr19 --output=chr19_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr19"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chr19 --output=chr19_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chr19"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chr19_average_ESC_100bp.txt --input2=chr19_average_MEF_100bp.txt --output1=chr19_MEF_less_ESC_0.99.txt --output2=chr19_MEF_more_ESC_0.99.txt --chromosome=chr19 --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chrX --output=chrX_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chrX"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chrX --output=chrX_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chrX"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chrX_average_ESC_100bp.txt --input2=chrX_average_MEF_100bp.txt --output1=chrX_MEF_less_ESC_0.99.txt --output2=chrX_MEF_more_ESC_0.99.txt --chromosome=chrX  --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3

perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_ESC/chrY --output=chrY_average_ESC_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chrY"
perl /storage/projects/nuctools.2.0/stable_nucs_replicates.pl --fileExtention=occ --inputDir=/storage/projects/mESC/100bp_MEFv2/chrY --output=chrY_average_MEF_100bp.txt --coordsCol=0 --occupCol=1 --StableThreshold=0.5 --chromosome="chrY"
perl /storage/projects/nuctools.2.0/compare_two_conditions.pl --input1=chrY_average_ESC_100bp.txt --input2=chrY_average_MEF_100bp.txt --output1=chrY_MEF_less_ESC_0.99.txt --output2=chrY_MEF_more_ESC_0.99.txt --chromosome=chrY  --windowSize=100 --threshold1=0.99 --threshold2=-0.99 --Col_coord=1 --Col_signal=3