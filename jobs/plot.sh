#!/bin/bash
#SBATCH --job-name=plot                         # Job name
#SBATCH --mail-type=NONE                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=heikki.tolppanen@oulu.fi    # Where to send mail	
#SBATCH --nodes=1                               # Run all processes on a single node	
#SBATCH --ntasks=1                              # Run a single task		
#SBATCH --cpus-per-task=4                       # Number of CPU cores per task
#SBATCH --mem=2gb                               # Job memory request
#SBATCH --time=0:10:00                          # Time limit hrs:min:sec
#SBATCH --output=batch_log_files/%x_%j_log.txt.txt     # Standard output and error log
pwd; hostname; date

/homedir02/htolppan/.julia/juliaup/julia-1.10.5+0.x64.linux.gnu/bin/julia -t 40 -e 'include("LehmusThreads.jl");
include("Physics/LeakageRemovalUnit/PlotsForPaper/3_figure/size_realistic/paper_plots/paper_plot.jl")'

date

# for 40 paral we would expect roughly 20h to calculate?
