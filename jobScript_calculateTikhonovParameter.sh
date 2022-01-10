#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=10.1G,s_vmem=10.0G

matlab -nodisplay -r "calculateTikhonovParameter; exit"
