# oxDNA run folder template

## Storage layout (same on DCC and Expanse; only the paths differ)
- **Home**: code, builds, conda envs, dotfiles only. No run folders.
- **SIM_ROOT** (linked as `~/sims`): run folders — inputs, job scripts, restart files,
  energy.dat, trajectory_cx1.dat, slurm logs. Work from `~/sims/<project>/<run>`.
  - DCC: `/hpc/dctrl/dduke/sims`
  - Expanse: `/expanse/lustre/projects/csd913/dduke/sims`
- **BULK_ROOT**: only the raw trajectory.dat (and its .bak copies).
  - DCC: `/work/dduke/bulk` (purged after 75 days)
  - Expanse: `/expanse/lustre/projects/csd913/dduke/bulk`
- Each run folder has a `scratch` link to its own bulk folder, made by run_oxdna.sh the
  first time it runs there. oxdna.in uses `trajectory_file = scratch/trajectory.dat`.
  Bulk folders are named `<slurm job ID>_<run folder name>` (e.g. `41234567_prod`), and
  each has a `run_folder.txt` saying which run made it.

## Files
- `job.sh` — Expanse job script. `job_dcc.sh` — DCC job script (rename to job.sh on DCC).
- `run_stages.sh` — runs each stage (min, relax, prod) in its own folder.
- `run_oxdna.sh` — runs one oxDNA stage; identical on both clusters.

## One-time setup (each cluster)
```bash
ln -s $SIM_ROOT ~/sims
mkdir -p $BULK_ROOT   # run_oxdna.sh expects it to exist
```
In `~/.bash_profile` (remove the old SCRATCH, MY_HOME, MY_HOME_LEN, HOME_LEN):
```bash
export BULK_ROOT=/work/dduke/bulk                                # DCC
export BULK_ROOT=/expanse/lustre/projects/csd913/dduke/bulk      # Expanse
rmsim() { rm -rf "$(readlink "$1/scratch")" "$1"; }   # delete a run and its bulk folder
```

## Usage
Pull to Mac (everything except the raw trajectory):
```bash
rsync -av --exclude='scratch' --exclude='*.tmp' dcc:sims/<project>/ ~/sims/<project>/
```
(or `expanse:` instead of `dcc:`)

Copy a run folder to start a new run without sharing its trajectory:
```bash
rsync -a --exclude=scratch old/ new/
```

Delete a run and its bulk folder:
```bash
rmsim <run_folder>
```

Find orphaned bulk folders (nothing links to them anymore):
```bash
find ~/sims/ -name scratch -type l -exec readlink {} \; | sort > linked.txt
ls -d $BULK_ROOT/* | sort > all.txt
comm -13 linked.txt all.txt
```
