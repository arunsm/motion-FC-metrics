#!/bin/bash

TASK=2
PATH=$PATH:$CFNAPPS/matlab/R2017a/bin/

ARRAY=( $(tail -n +2 /data/jag/bassett-lab/hcp_Max/Data/Covariates/S1200_Release_Subjects_Demographics.csv | cut -d , -f1))

LENGTH=${#ARRAY[@]}
echo Total number of subjects: $LENGTH

#$ -t 1-1206
#$ -j y
#$ -cwd
#$ -q himem.q,all.q,basic.q,gpu.q
#$ -l h_vmem=6.1G,s_vmem=6.0G

i=`expr ${SGE_TASK_ID} - 1`

if [[ $i -ge ${LENGTH} ]]; then
  echo 'Array index > than number of elements'
else
  SUB=${ARRAY[$i]}
  echo Submitting job for subject $SUB
  ./jobScript_FC_perSubject.sh $SUB $TASK 
fi
