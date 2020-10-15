#!/usr/bin/env bash

task=rfMRI_REST1_LR;
mkdir /Users/arunmahadevan/Documents/HCP_motionData/$task
ARRAY=( $(tail -n +2 /Users/arunmahadevan/Dropbox/Mahadevan_Bassett_Projects/Motion_FunctionalConnectivity/data/Covariates/S1200_Release_Subjects_Demographics.csv | cut -d , -f1))

for SUB in ${ARRAY[@]}
	do
	cp /Users/arunmahadevan/.CMVolumes/asm7/hcp-openaccess/HCP_1200/$SUB/MNINonLinear/Results/$task/Movement_Regressors.txt /Users/arunmahadevan/Documents/HCP_motionData/$task
	mv /Users/arunmahadevan/Documents/HCP_motionData/$task/Movement_Regressors.txt /Users/arunmahadevan/Documents/HCP_motionData/$task/${SUB}_Movement_Regressors.txt
	echo $SUB
	done
