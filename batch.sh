#!/bin/bash
#SBATCH --job-name=fig2                         # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=heikki.tolppanen@oulu.fi    # Where to send mail	
#SBATCH --nodes=1                               # Run all processes on a single node	
#SBATCH --ntasks=1                              # Run a single task		
#SBATCH --cpus-per-task=40                      # Number of CPU cores per task
#SBATCH --mem=12gb                              # Job memory request
#SBATCH --time=18:00:00                         # Time limit hrs:min:sec
#SBATCH --output=batch_run.log                  # Standard output and error log
pwd; hostname; date

/homedir02/htolppan/.julia/juliaup/julia-1.10.5+0.x64.linux.gnu/bin/julia -t 40 -e 'include("LehmusJulia.jl");
include("Physics/LeakageRemovalUnit/PlotsForPaper/second_figure/total_anh/gamma/code.jl")'

date

# for 40 paral we would expect roughly 20h to calculate?
