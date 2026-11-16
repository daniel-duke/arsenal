#!/bin/bash

# Template version: 2.0 on 11/11/26

## Description
# this script runs oxDNA with the input file in the current folder, using the
  # GPU or CPU build of oxDNA depending on the backend set by the input file.
# the raw trajectory is written to scratch/trajectory.dat and periodically
  # coarsened to the run folder.
# bulk trajectories are stored in the scratch folder; if this folder (or link)
  # does not exist, it is created as follows: if SCRATCH_ROOT is set, scratch
  # is a link to a new, randomly named folder in SCRATCH_ROOT; otherwise, scratch
  # is a plain folder in the run folder.

### environment
oxdna_gpu=$OXDNA_GPU
oxdna_cpu=${OXDNA_CPU:-$OXDNA_GPU}
coarsener=$ars/coarsenOxDNA.py
scratch_root=$SCRATCH_ROOT

### settings
inFile=oxdna.in
coarsen_args="--coarse-time 1 --precision 3"
coarsen_every=60 # seconds

### kick things off
echo "Start Time: $(date)"
echo -e "Directory: $(pwd)\n"
START_EPOCH=$(date +%s)
set -e

### get parameter value from input file
getParam() {
	sed -n "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*\([^[:space:]#]*\).*/\1/p" "$1" | tail -n 1
}

### choose oxdna executable from backend
backend=$(getParam $inFile backend)
if [ "$backend" = "CUDA" ]; then
	oxdna=$oxdna_gpu
elif [ "$backend" = "CPU" ]; then
	oxdna=$oxdna_cpu
fi
echo "oxDNA: $oxdna"

### set up scratch folder
if [[ ! -d scratch ]]; then
	if [[ -n $scratch_root ]]; then
		scratch_target=$(mktemp -d "$scratch_root/XXXXXXXXXX")
		echo "$PWD" > "$scratch_target/run_folder.txt"
		ln -sfn "$scratch_target" scratch
	else
		mkdir scratch
	fi
fi
echo "Scratch: $(readlink scratch || echo "local folder")"

### record input
cat $inFile

### run oxDNA in background
$oxdna $inFile &
oxdna_pid=$!

### wait and coarsen loop
while kill -0 $oxdna_pid 2>/dev/null; do

	# wait for coarsen interval, checking every second if oxDNA has finished
	for ((t=0; t<coarsen_every; t+=1)); do
		kill -0 $oxdna_pid 2>/dev/null || break
		sleep 1
	done
	kill -0 $oxdna_pid 2>/dev/null || break

	# coarsen silently, ignoring errors
	python $coarsener scratch/trajectory.dat $coarsen_args --outHere --update --hush > /dev/null 2>&1 || true
done

### get oxdna exit status
oxdna_status=0
wait $oxdna_pid || oxdna_status=$?
echo ""

### final coarsening, heeding errors
python $coarsener scratch/trajectory.dat $coarsen_args --outHere --update || true

### wrap things up
END_EPOCH=$(date +%s)
echo "Job completed: $(date) - $(($END_EPOCH - $START_EPOCH)) seconds elapsed"
exit $oxdna_status
