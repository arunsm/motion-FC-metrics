#!/bin/bash
#$ -j y
#$ -l h_vmem=32.0G,s_vmem=31.0G

IMAGE_DIR=$HOME/singularityImages
VERSION="latest"

singularity build $IMAGE_DIR/octave-$VERSION.simg docker://mtmiller/octave:$VERSION
