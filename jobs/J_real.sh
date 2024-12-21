#!/bin/bash
#SBATCH --job-name=J_real                       # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=heikki.tolppanen@oulu.fi    # Where to send mail	
#SBATCH --nodes=1                               # Run all processes on a single node	
#SBATCH --ntasks=1                              # Run a single task		
#SBATCH --cpus-per-task=40                      # Number of CPU cores per task
#SBATCH --mem=20gb                              # Job memory request
#SBATCH --time=3:00:00                         # Time limit hrs:min:sec
#SBATCH --output=jobs/logs/%x_%j_log.txt     # Standard output and error log

file_to_run="Physics/LeakageRemovalUnit/PlotsForPaper/3_figure/J_realistic/make_T_comparisons.jl"

echo "$file_to_run"; hostname; date

julia -t 40 -e "
include(\"LehmusThreads.jl\");
include(\"$file_to_run\")
"

date