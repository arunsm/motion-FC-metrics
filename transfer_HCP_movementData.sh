#!/bin/bash/

task=rfMRI_REST2_RL;
mkdir /Users/ArunMahadevan/Documents/hcp_Max/Motion_S1200/$task
ARRAY=( $(tail -n +2 /Users/ArunMahadevan/Documents/hcp_Max/Covariates/S1200_Release_Subjects_Demographics.csv | cut -d , -f1))

for SUB in ${ARRAY[@]}
	do
	cp /Users/ArunMahadevan/.CMVolumes/asm7/HCP_1200/$SUB/MNINonLinear/Results/$task/Movement_RelativeRMS_mean.txt /Users/ArunMahadevan/Documents/hcp_Max/Motion_S1200/$task
	mv /Users/ArunMahadevan/Documents/hcp_Max/Motion_S1200/$task/Movement_RelativeRMS_mean.txt /Users/ArunMahadevan/Documents/hcp_Max/Motion_S1200/$task/${SUB}_Movement_RelativeRMS_mean.txt
	echo $SUB
	done
