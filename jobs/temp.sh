#!/bin/bash
#SBATCH --job-name=fig2_dt_tes                  # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=heikki.tolppanen@oulu.fi    # Where to send mail	
#SBATCH --nodes=1                               # Run all processes on a single node	
#SBATCH --ntasks=1                              # Run a single task		
#SBATCH --cpus-per-task=40                      # Number of CPU cores per task
#SBATCH --mem=60gb                              # Job memory request
#SBATCH --time=30:00:00                         # Time limit hrs:min:sec
#SBATCH --output=jobs/logs/%x_%j_log.txt        # Standard output and error log

file_to_run="Physics/LeakageRemovalUnit/PlotsForPaper/3_figure/size/size_dependance_optimal_spot/code.jl"
traj="2000"

echo "$file_to_run"; echo "trajectories = $traj"; hostname; date; echo

start_time=$(date +%s)  # Record start time
julia -t 40 -e "
ARGS = [\"$traj\"]
include(\"LehmusThreads.jl\");
include(\"$file_to_run\")
"
end_time=$(date +%s)    # Record end time
elapsed_time=$((end_time - start_time))  # Calculate elapsed time

echo
echo "Trajectories: $traj"
./jobs/_elapsed_time.sh $elapsed_time
date
