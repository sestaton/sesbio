

qstat -u "*" | grep "mri-q32-30d" | awk -F' ' '{print $NF}' | awk '{sum+=$1} END {print sum}'