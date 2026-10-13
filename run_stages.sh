# Template version: 1.0 on 10/09/26

## Description
# this script runs a series of oxdna simulation stages, each in its own folder.
# for each stage, it creates a folder named after the stage, copies in the
  # topology, the stage input file ({stage}.in, renamed to oxdna.in),
  # run_oxdna.sh, and the starting configuration, then runs the stage and
  # waits for it to finish before moving on to the next one.
# the first stage starts from the initial configuration file; every later
  # stage starts from the final configuration of the previous stage, which is
  # copied in as last_conf_{previous stage}.dat.
# the conf_file in each stage input file must match its starting configuration
  # name, and the topology must match the topology file name (both checked at
  # the start).

stages=(min relax)
confFile=square_v3-5_ideal.dat
topFile=square_v3-5.top
lastConfFile=last_conf.dat  # final configuration written by oxdna

set -e
echo "Stages started: $(date)"
echo "Directory: $(pwd)"
START_EPOCH=$(date +%s)

### get value of parameter from oxdna input file (empty if not set)
getParam() {
    sed -n "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*\([^[:space:]#]*\).*/\1/p" "$1" | tail -n 1
}

### check for required files
error=0
for file in "$confFile" "$topFile" run_oxdna.sh "${stages[@]/%/.in}"; do
    if [ ! -f "$file" ]; then
        echo "Error: Could not find required file: $file"
        error=1
    fi
done
[ $error -eq 0 ] || exit 1

### check stage folders are not already there
for stage in "${stages[@]}"; do
    if [ -e "$stage" ]; then
        echo "Error: Stage folder already exists: $stage"
        error=1
    fi
done

### check input files point to correct files
startConf=$confFile
for stage in "${stages[@]}"; do
    inFile=$stage.in
    if [ "$(getParam $inFile conf_file)" != "$startConf" ]; then
        echo "Error: conf_file in $inFile must be $startConf."
        error=1
    fi
    if [ "$(getParam $inFile topology)" != "$topFile" ]; then
        echo "Error: topology in $inFile must be $topFile."
        error=1
    fi
    lastConf=$(getParam $inFile lastconf_file)
    if [ -n "$lastConf" ] && [ "$lastConf" != "$lastConfFile" ]; then
        echo "Error: lastconf_file in $inFile must be $lastConfFile (or unset)."
        error=1
    fi
    startConf=last_conf_$stage.dat
done
[ $error -eq 0 ] || exit 1

### run stages
prevStage=""
for stage in "${stages[@]}"; do

    ### populate stage folder
    mkdir "$stage"
    cp "$topFile" run_oxdna.sh "$stage/"
    cp "$stage.in" "$stage/oxdna.in"
    if [ -z "$prevStage" ]; then
        cp "$confFile" "$stage/"
    else
        cp "$prevStage/$lastConfFile" "$stage/last_conf_$prevStage.dat"
    fi

    ### run stage and wait for it to finish
    echo ""
    echo "Stage $stage started: $(date)"
    status=0
    (cd "$stage" && bash run_oxdna.sh >> report.out 2>&1) || status=$?
    if [ $status -ne 0 ]; then
        echo "Error: Stage $stage failed with status $status (see $stage/report.out)."
        exit 1
    fi
    if [ ! -f "$stage/$lastConfFile" ]; then
        echo "Error: Stage $stage did not write $lastConfFile."
        exit 1
    fi
    echo "Stage $stage completed: $(date)"
    prevStage=$stage
done

END_EPOCH=$(date +%s)
echo ""
echo "Stages completed: $(date) - $(($END_EPOCH - $START_EPOCH)) seconds elapsed"