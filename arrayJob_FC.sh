#!/bin/bash
#$ -cwd

TASK=4

ARRAY=( $(tail -n +2 /cbica/home/mahadeva/motion-FC-metrics/data/Covariates/S1200_Release_Subjects_Demographics.csv | cut -d , -f1))

LENGTH=${#ARRAY[@]}
echo Total number of subjects: $LENGTH

#$ -t 1-1206
#$ -j y
#$ -l h_vmem=6.1G,s_vmem=6.0G

i=`expr ${SGE_TASK_ID} - 1`

if [[ $i -ge ${LENGTH} ]]; then
  echo 'Array index > than number of elements'
else
  SUB=${ARRAY[$i]}
  echo Submitting job for subject $SUB
  ./jobScript_FC_perSubject_gordon.sh $SUB $TASK
  ./jobScript_FC_perSubject_yeo.sh $SUB $TASK 
fi
