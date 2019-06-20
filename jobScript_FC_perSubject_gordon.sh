#!/bin/bash

subjectID=$1
task=$2
echo $task
matlab -nodisplay -r "engine_FunctionalConnectivity_gordon($subjectID, $task); exit"
