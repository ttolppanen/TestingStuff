#!/bin/bash
#SBATCH --job-name=J_real                       # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=heikki.tolppanen@oulu.fi    # Where to send mail	
#SBATCH --nodes=1                               # Run all processes on a single node	
#SBATCH --ntasks=1                              # Run a single task		
#SBATCH --cpus-per-task=20                      # Number of CPU cores per task
#SBATCH --mem=80gb                              # Job memory request
#SBATCH --time=50:00:00                         # Time limit hrs:min:sec
#SBATCH --output=jobs/logs/%x_%j_log.txt     # Standard output and error log

file_to_run="Physics/LeakageRemovalUnit/PlotsForPaper/appendix/3_figure_dissipation_high/J_realistic/make_T_comparisons.jl"
traj="1000" # 1000 traj 45min -> 1333/h -> 40k takes 30h!

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

