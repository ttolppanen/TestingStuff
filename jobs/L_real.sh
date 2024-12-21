#!/bin/bash
#SBATCH --job-name=L_real_T_times               # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=heikki.tolppanen@oulu.fi    # Where to send mail	
#SBATCH --nodes=1                               # Run all processes on a single node	
#SBATCH --ntasks=1                              # Run a single task		
#SBATCH --cpus-per-task=40                      # Number of CPU cores per task
#SBATCH --mem=40gb                              # Job memory request
#SBATCH --time=12:00:00                         # Time limit hrs:min:sec
#SBATCH --output=jobs/logs/%x_%j_log.txt     # Standard output and error log
pwd; hostname; date

/homedir02/htolppan/.julia/juliaup/julia-1.10.5+0.x64.linux.gnu/bin/julia -t 40 -e 'include("LehmusThreads.jl");
include("Physics/LeakageRemovalUnit/PlotsForPaper/3_figure/size_realistic/make_T_comparisons.jl")'

date
