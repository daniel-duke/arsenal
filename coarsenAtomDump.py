import numpy as np
import argparse
import sys
import os
import re
import shutil

## Description
# this script reads a LAMMPS-style trajectory file (or any dump file with the
  # appropriate header) and rewrites the data into a new file with coarsened
  # time steps.
# this script runs from the terminal with the trajectory file to coarsen as its
  # argument; the coarsening is set by one of two options: (1) coarse-time
  # keeps every nth frame, regardless of step number, and (2) steps-per-frame
  # keeps only frames whose step number is divisible by it.
# step numbers are kept from the original trajectory, unless --setSteps is
  # given (only with --steps-per-frame), which renumbers the starting from 0.
# dryRun reports the step counts without writing anything; adding printSteps
  # also lists the step of every frame and its step after coarsening.


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('datFile',type=str)
	parser.add_argument('--out',type=str,default=None)
	parser.add_argument('--coarse-time',type=int,default=None)
	parser.add_argument('--steps-per-frame',type=int,default=None)
	parser.add_argument('--setSteps',action='store_true')
	parser.add_argument('--update',action='store_true')
	parser.add_argument('--dryRun',action='store_true')
	parser.add_argument('--printSteps',action='store_true')
	parser.add_argument('--hush',action='store_true')
	args = parser.parse_args()

	### check arguments
	if args.coarse_time is not None and args.steps_per_frame is not None:
		print("Error: Cannot set both coarse-time and steps-per-frame.\n")
		sys.exit(1)
	if args.setSteps and args.steps_per_frame is None:
		print("Error: setSteps can only be used with steps-per-frame.\n")
		sys.exit(1)
	if args.printSteps and not args.dryRun:
		print("Error: printSteps can only be used with dryRun.\n")
		sys.exit(1)
	if args.coarse_time is not None and args.coarse_time < 1:
		print("Error: Coarse time must be positive.\n")
		sys.exit(1)
	if args.steps_per_frame is not None and args.steps_per_frame < 1:
		print("Error: Steps per frame must be positive.\n")
		sys.exit(1)
	if args.coarse_time is None and args.steps_per_frame is None:
		args.coarse_time = 1

	### read old, write new
	readWriteAtomDump(args.datFile, args.out, args.coarse_time, args.steps_per_frame, args.setSteps, args.hush, args.update, args.dryRun, args.printSteps)
	if not args.hush: print()


################################################################################
### File Managers

### read LAMMPS-style atom dump, write coarsened one
def readWriteAtomDump(datFile, outDatFile, coarse_time, steps_per_frame, setSteps, hush, update, dryRun, printSteps):
	if outDatFile is None:
		outDatFile = getCoarseName(datFile, coarse_time, steps_per_frame)

	### extract metadata
	checkFileExist(datFile, "trajectory", requireData=True)
	with open(datFile, 'r') as f:

		for i in range(2): line = f.readline()
		step_initial = int(line.split()[0])

		for i in range(2): line = f.readline()
		nbd_total = int(line.split()[0])

		for i in range(nbd_total+7): line = f.readline()
		if not line:
			print("Error: Trajectory needs at least two frames.\n")
			sys.exit(1)
		steps_per_frame_initial = int(line.split()[0]) - step_initial

	### check step compatibility
	if steps_per_frame is not None and steps_per_frame % steps_per_frame_initial != 0:
		print(f"Flag: Steps per frame is not a multiple of the initial steps per frame ({steps_per_frame_initial}).")

	### report step counts
	if not hush:
		with open(datFile, 'rb') as f:
			nstep_recorded = sum(1 for line in f if line.startswith(b'ITEM: TIMESTEP'))
		nstep_coarse = countCoarseFrames(step_initial, steps_per_frame_initial, nstep_recorded, coarse_time, steps_per_frame, 0)
		nstep_simulation = (nstep_recorded-1)*steps_per_frame_initial
		print("{:1.2e} steps in simulation".format(nstep_simulation))
		print("{:1.2e} steps in trajectory".format(nstep_recorded))
		print("{:1.2e} steps after coarsening".format(nstep_coarse))

	### open source trajectory
	with open(datFile, 'r') as fin:
		frames = coarsenFrames(readFrames(fin), coarse_time, steps_per_frame, setSteps, 0)

		### check for old coarse file
		nstep_old = 0
		if update and os.path.isfile(outDatFile):
			nstep_old = checkOldFile(outDatFile, frames)
			if not hush:
				print("{:1.2e} steps already coarsened".format(nstep_old))
		else:
			update = False

		### stop before changing anything
		if dryRun:
			if printSteps:
				printFrameSteps(datFile, coarse_time, steps_per_frame, setSteps)
			return

		### write coarse trajectory to temporary file
		tmpDatFile = outDatFile + ".tmp"
		with open(tmpDatFile, 'w', buffering=1024*1024) as fout:

			### copy the existing coarse frames
			if update:
				with open(outDatFile, 'r') as fold:
					shutil.copyfileobj(fold, fout)

			### then add the new ones
			for i, frame_lines in enumerate(frames, nstep_old):

				### keep in touch (status bar only appears if there are new steps)
				if not hush:
					if i == nstep_old:
						initStatusBar("Coarsening trajectory")
						for j in range(nstep_old):
							updateStatusBar(j, nstep_coarse)
					updateStatusBar(i, nstep_coarse)

				### copy frame
				fout.write(''.join(frame_lines))

		### swap in the finished file
		os.replace(tmpDatFile, outDatFile)


### check that old coarsened file matches the start of the coarsened trajectory
def checkOldFile(outDatFile, frames):
	nstep_old = 0

	### open old file
	with open(outDatFile, 'r') as fold:

		### loop over old frames
		for frame_old in readFrames(fold):

			### get matching new frame
			frame_new = next(frames, None)
			if frame_new is None:
				print(f"Error: Old file has more time steps than new file.\n")
				sys.exit(1)

			### compare timesteps
			if frame_old[1] != frame_new[1]:
				print(f"Error: Time step mismatch between old and new file.\n")
				sys.exit(1)

			### compare first entry (if any)
			if frame_old[9:10] != frame_new[9:10]:
				print(f"Error: Atom data mismatch between old and new file.\n")
				sys.exit(1)

			### update counter
			nstep_old += 1

	### results
	return nstep_old


### iterate over complete frames
def readFrames(fin):
	while True:

		### check for new frame
		line = fin.readline()
		if not line:
			return

		### read header through entry count
		frame_lines = [line] + [fin.readline() for _ in range(3)]
		if not frame_lines[-1].endswith('\n'):
			return

		### read rest of frame (entry count can change between frames)
		nentry = int(frame_lines[3].split()[0])
		frame_lines += [fin.readline() for _ in range(nentry+5)]

		### stop at incomplete frame (simulation still writing)
		if not frame_lines[-1].endswith('\n'):
			return

		yield frame_lines


### print the step of every frame, next to its step after coarsening
def printFrameSteps(datFile, coarse_time, steps_per_frame, setSteps):
	print(f"{'original':<12} coarsened")
	with open(datFile, 'r') as fin:
		for frame_lines, step, step_new in labelFrames(readFrames(fin), coarse_time, steps_per_frame, setSteps, 0):
			print(f"{step:<12} {'' if step_new is None else step_new}")


### keep the frames that survive coarsening, optionally renumbering their steps
def coarsenFrames(frames, coarse_time, steps_per_frame, setSteps, offset):
	for frame_lines, step, step_new in labelFrames(frames, coarse_time, steps_per_frame, setSteps, offset):
		if step_new is None:
			continue
		if step_new != step:
			frame_lines[1] = f"{step_new}\n"
		yield frame_lines


### pair each frame with its step and its step after coarsening (None if dropped)
def labelFrames(frames, coarse_time, steps_per_frame, setSteps, offset):
	nkept = 0
	for i, frame_lines in enumerate(frames):
		step = int(frame_lines[1].split()[0])

		### mark unused frames
		if coarse_time is not None:
			keep = (i+offset) % coarse_time == 0
		else:
			keep = step % steps_per_frame == 0
		if not keep:
			yield frame_lines, step, None
			continue

		### renumber steps
		step_new = (nkept+offset)*steps_per_frame if setSteps else step
		nkept += 1
		yield frame_lines, step, step_new


################################################################################
### Utility Functions

### determine if file exists
def checkFileExist(file, name="the", required=True, requireData=False, hush=False):
	if os.path.isfile(file):
		if not requireData or os.path.getsize(file):
			return True
		else:
			print(f"Error: {name} file has no content:")
			print(file + "\n")
			sys.exit(1)
	else:
		if required:
			print(f"Error: Could not find {name} file:")
			print(file + "\n")
			sys.exit(1)
		else:
			if not hush:
				print(f"Flag: Could not find {name} file.")
			return False


### count frames that survive coarsening (offset is 1 if step 0 is not recorded)
def countCoarseFrames(step_initial, steps_per_frame_initial, nstep_recorded, coarse_time, steps_per_frame, offset):
	if coarse_time is not None:
		return (nstep_recorded-1+offset) // coarse_time - (offset-1) // coarse_time
	return sum(1 for i in range(nstep_recorded) if (step_initial+i*steps_per_frame_initial) % steps_per_frame == 0)


### name coarse file, combining consecutive coarse-time factors
def getCoarseName(file, coarse_time, steps_per_frame):
	base, ext = os.path.splitext(file)
	if steps_per_frame is not None:
		return f"{base}_spf{steps_per_frame}{ext}"
	_, sep, suffix = os.path.basename(base).rpartition('_')
	match = re.fullmatch(r'cx(\d+)', suffix) if sep else None
	if match:
		return f"{base[:-len(suffix)]}cx{int(match.group(1))*coarse_time}{ext}"
	return f"{base}_cx{coarse_time}{ext}"


### print status bar header
def initStatusBar(description, length=30):

	### dynamic output
	if sys.stdout.isatty():
		print(f"{description}:")

	### static output
	else:
		if len(description) > length-8:
			print("Flag: Status bar description too long.")
			description = description[:length-6]
		pad = length - len(description) - 4
		print(f"{'='*(length)}")
		print(f"-- {description} {'-'*pad}")


### print status bar during loop iterations
def updateStatusBar(i, n, length=30, units='steps'):

	### dynamic output
	if sys.stdout.isatty():
		if i < n-1:
			print(f"\r-- {i+1}/{n} {units}", end='', flush=True)
		else:
			print(f"\r-- {n}/{n} {units}")

	### static output
	else:
		nchar_before = length*i // n
		nchar_after = length*(i+1) // n
		nchar_add = nchar_after - nchar_before
		print("=" * nchar_add, end='', flush=True)
		if i == n-1:
			print(flush=True)


### run the script
if __name__ == "__main__":
	main()

