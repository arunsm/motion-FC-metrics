#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=6.1G,s_vmem=6.0G

matlab -nodisplay -r "computeQCFC; exit"
