set -euo pipefail
#$ -j y
#$ -o /cbica/home/tooleyu/projects/in_progress/arun_fc_metrics_motion/
#$ -l h_vmem=10.5G,s_vmem=10.3G
#$ -V
#$ -cwd

code_dir='/cbica/home/tooleyu/projects/in_progress/motion-FC-metrics/system_identifiability_analyses'
matlab -nodisplay -r "cd ${code_dir}, run('modularity_quality_analysis.m'); exit"
