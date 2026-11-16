#!/bin/bash

# Template version: 2.0 on 11/11/26

## Description
# this script runs a series of oxdna simulation stages, each in its own folder,
  # all within one slurm job.
# for each stage, it creates a folder named after the stage, copies in the
  # appropriate files (topology, configuration, stage input file ({stage}.in,
  # renamed to oxdna.in), and run_oxdna.sh), and runs the stage.
# the topology, initial configuration, and final configuration parameters are
  # set automatically in each stage's input file, with a warning printed if the
  # input file has any value other than the correct one or BLANK.
# the first stage starts from the initial configuration, and every later
  # stage starts from the final configuration of the previous stage, which is
  # copied in as last_conf_{previous stage}.dat.
# if a file called forces.txt exists, it is copied into every stage folder;
  # otherwise, if forces_{stage}.txt exists, it is copied into that stage
  # folder as forces.txt.
# any existing stage folders are removed at the start, so a requeued job
  # (for example, after preemption) starts over from the first stage.

### parameters
stages=(min relax prod)
topFile=ds_bp256.top
confFile=ds_bp256_ideal.dat

### kick things off
echo "Start Time: $(date)"
echo -e "Directory: $(pwd)\n"
START_EPOCH=$(date +%s)
set -e

### get parameter from input file
getParam() {
	sed -n "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*\([^[:space:]#]*\).*/\1/p" "$1" | tail -n 1
}

### check for correct parameter in input file
checkParam() {
	local current=$(getParam "$1" "$2")
	if [ -z "$current" ]; then
		echo "Error: $2 not set in $1."
		error=1
	elif [ "$current" != "$3" ] && [ "$current" != "BLANK" ]; then
		echo "Warning: $2 in $1 is $current, will replace it with $3."
	fi
}

### set parameter in input file
setParam() {
	sed "s|^\([[:space:]]*$2[[:space:]]*=[[:space:]]*\)[^[:space:]#]*|\1$3|" "$1" > "$1.tmp"
	mv "$1.tmp" "$1"
}

### check current folder for required files
error=0
for file in "$confFile" "$topFile" run_oxdna.sh "${stages[@]/%/.in}"; do
	if [ ! -f "$file" ]; then
		echo "Error: Could not find required file: $file"
		error=1
	fi
done
[ $error -eq 0 ] || exit 1

### check input file for correct file references
confFileStart=$confFile
for stage in "${stages[@]}"; do
	checkParam "$stage.in" topology "$topFile"
	checkParam "$stage.in" conf_file "$confFileStart"
	checkParam "$stage.in" lastconf_file last_conf.dat
	confFileStart=last_conf_$stage.dat
done
[ $error -eq 0 ] || exit 1

### remove existing stage folders and their scratch folders
for stage in "${stages[@]}"; do
	if [ -e "$stage" ]; then
		echo "Warning: Stage folder '$stage' already exists, removing it and its scratch folder."
		rm -rf "$(readlink "$stage/scratch")" "$stage"
	fi
done

### run stages
stage_prev=""
for stage in "${stages[@]}"; do

	### populate stage folder
	mkdir "$stage"
	cp "$topFile" run_oxdna.sh "$stage/"
	cp "$stage.in" "$stage/oxdna.in"
	if [ -f forces.txt ]; then
		cp forces.txt "$stage/forces.txt"
	elif [ -f "forces_$stage.txt" ]; then
		cp "forces_$stage.txt" "$stage/forces.txt"
	fi
	if [ -z "$stage_prev" ]; then
		confFileStart=$confFile
		cp "$confFile" "$stage/"
	else
		confFileStart=last_conf_$stage_prev.dat
		cp "$stage_prev/last_conf.dat" "$stage/$confFileStart"
	fi

	### set file parameters in input file
	setParam "$stage/oxdna.in" topology "$topFile"
	setParam "$stage/oxdna.in" conf_file "$confFileStart"
	setParam "$stage/oxdna.in" lastconf_file last_conf.dat

	### run stage and wait for it to finish
	echo -e "Stage $stage started: $(date)"
	status=0
	(cd "$stage" && bash run_oxdna.sh >> report.out 2>&1) || status=$?
	if [ $status -ne 0 ]; then
		echo "Error: Stage $stage failed with status $status (see $stage/report.out)."
		exit 1
	fi
	if [ ! -f "$stage/last_conf.dat" ]; then
		echo "Error: Stage $stage did not write last_conf.dat."
		exit 1
	fi
	echo -e "Stage $stage completed: $(date)\n"
	stage_prev=$stage
done

END_EPOCH=$(date +%s)
echo -e "Job completed: $(date) - $(($END_EPOCH - $START_EPOCH)) seconds elapsed"
