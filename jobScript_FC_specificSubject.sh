#!/bin/bash

# script to compute FC measures for a particular subject

#$ -cwd
#$ -l h_vmem=6.1G,s_vmem=6.0G

subjectID=580751
task=2
matlab -nodisplay -r "engine_FunctionalConnectivity_gordon($subjectID, $task); exit"
