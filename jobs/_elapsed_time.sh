#!/bin/bash
# _elapsed_time.sh: Calculate and format elapsed time

elapsed_time=$1  # Accept elapsed time in seconds as an argument

# Calculate hours, minutes, and seconds
hours=$((elapsed_time / 3600))
remainder_seconds=$((elapsed_time % 3600))
minutes=$((remainder_seconds / 60))
seconds=$((remainder_seconds % 60))

# Output the formatted elapsed time
echo "Elapsed time: ${hours}h ${minutes}min ${seconds}s"
