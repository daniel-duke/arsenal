#!/bin/bash

# Template version: 3.0 on 11/19/26

## Description
# this script runs LAMMPS with the input file in the current folder, using MPI
  # tasks and/or OpenMP threads as set below (0 turns either off).
# the raw trajectory is written to trajectory.dat and periodically centered
  # to trajectory_centered.dat.

### environment
lammps=$LAMMPS
centerer=$ARS/centerAtomDump.py

### settings
inFile=lammps.in
nthread=0
ntask=1
center_args="--center com --unwrap"
center_every=60 # seconds

### kick things off
echo "Start Time: $(date)"
echo -e "Directory: $(pwd)\n"
START_EPOCH=$(date +%s)
set -e

### check lammps executable
if [ -z "$lammps" ] || [ ! -f "$lammps" ]; then
	echo "Error: Could not find LAMMPS executable."
	exit 1
fi
echo "LAMMPS: $lammps"

### build run command from parallel settings
cmd=("$lammps" -in $inFile -log none)
if [ $nthread -gt 0 ]; then
	cmd+=(-sf omp -pk omp $nthread)
fi
if [ $ntask -gt 0 ]; then
	cmd=(mpirun -n $ntask "${cmd[@]}")
fi

### record input
cat $inFile

### run LAMMPS in background
"${cmd[@]}" &
lammps_pid=$!

### wait and center loop
while kill -0 $lammps_pid 2>/dev/null; do

	# wait for center interval, checking every second if LAMMPS has finished
	for ((t=0; t<center_every; t+=1)); do
		kill -0 $lammps_pid 2>/dev/null || break
		sleep 1
	done
	kill -0 $lammps_pid 2>/dev/null || break

	# center silently, ignoring errors
	python $centerer trajectory.dat $center_args --hush > /dev/null 2>&1 || true
done

### get lammps exit status
lammps_status=0
wait $lammps_pid || lammps_status=$?
echo ""

### final centering, heeding errors
python $centerer trajectory.dat $center_args || true

### wrap things up
END_EPOCH=$(date +%s)
echo "Job completed: $(date) - $(($END_EPOCH - $START_EPOCH)) seconds elapsed"
exit $lammps_status
