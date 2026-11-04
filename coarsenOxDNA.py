import numpy as np
import argparse
import sys
import os
import re
import shutil

## Description
# this script reads an oxDNA trajectory file and rewrites the data into a new
  # file with coarsened time steps.
# this can help reduce file sizes when the output frequency of the simulation
  # was not set appropriately, or when the simulation needed to be run for
  # longer than expected.
# this script loads the old file and writes the new file line by line, which 
  # can be slower for small scripts, but it never times out.
# this script runs from the terminal with the trajectory file to coarsen as its
  # argument; the coarsening is set by one of two options:
  # --coarse-time keeps every nth frame, regardless of step number (oxDNA does
  # not record step 0, so the kept steps are n, 2n, ... times the original
  # steps per frame, just as if oxDNA had been run with the coarser output).
  # --steps-per-frame keeps only frames whose step number is divisible by it.
# step numbers are kept from the original trajectory, unless --setSteps is
  # given (only with --steps-per-frame), which renumbers the kept frames as
  # steps_per_frame, 2*steps_per_frame, ...
# the output goes next to the input (with "_cx{coarse_time}" or
  # "_spf{steps_per_frame}" added to the name) unless --out is given;
  # coarsening an already coarse-timed file by coarse-time combines the
  # factors (e.g. "_cx2" then --coarse-time 3 gives "_cx6").
# --dryRun reports the step counts without writing anything; adding
  # --printSteps also lists the step of every frame and its step after
  # coarsening (blank if the frame is dropped).
# the output is written to a temporary file and then swapped in, so anything
  # reading the output (e.g. rsync) never sees a half-written file.


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('datFile',type=str)
	parser.add_argument('--out',type=str,default=None)
	parser.add_argument('--coarse-time',type=int,default=None)
	parser.add_argument('--steps-per-frame',type=int,default=None)
	parser.add_argument('--precision',type=int,default=8)
	parser.add_argument('--setSteps',action='store_true')
	parser.add_argument('--keepVelocity',action='store_true')
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
		print("Error: Coarse time must be 1 or greater.\n")
		sys.exit(1)
	if args.steps_per_frame is not None and args.steps_per_frame < 1:
		print("Error: Steps per frame must integer over 1.\n")
		sys.exit(1)
	if args.coarse_time is None and args.steps_per_frame is None:
		args.coarse_time = 1

	### read old, write new
	readWriteOxDNA(args.datFile, args.out, args.coarse_time, args.steps_per_frame, args.setSteps, args.hush, args.precision, args.keepVelocity, args.update, args.dryRun, args.printSteps)
	if not args.hush: print()


################################################################################
### File Managers

### read oxdna trajectory, write coarsened one
def readWriteOxDNA(datFile, outDatFile, coarse_time, steps_per_frame, setSteps, hush, precision, keepVelocity, update, dryRun, printSteps):
	if outDatFile is None:
		outDatFile = getCoarseName(datFile, coarse_time, steps_per_frame)

	### extract metadata
	checkFileExist(datFile, "trajectory", requireData=True)
	with open(datFile, 'r') as f:

		line = f.readline()
		step_initial = int(line.split()[2])

		line = f.readline()
		dbox = float(line.split()[2])

		line = f.readline()
		nba_total = 0
		line = f.readline()
		while line and line[0] != 't':
			nba_total += 1
			line = f.readline()
		if not line:
			print("Error: Trajectory needs at least two frames.\n")
			sys.exit(1)
		steps_per_frame_initial = int(line.split()[2]) - step_initial

	### check step compatibility
	if steps_per_frame is not None and steps_per_frame % steps_per_frame_initial != 0:
		print(f"Flag: Steps per frame is not a multiple of the initial steps per frame ({steps_per_frame_initial}).")

	### report step counts
	if not hush:
		with open(datFile, 'rb') as f:
			nline = sum(1 for _ in f)
		nstep_recorded = nline // (nba_total+3)
		nstep_coarse = countCoarseFrames(step_initial, steps_per_frame_initial, nstep_recorded, coarse_time, steps_per_frame, 1)
		nstep_simulation = nstep_recorded*steps_per_frame_initial
		print("{:1.2e} steps in simulation".format(nstep_simulation))
		print("{:1.2e} steps in trajectory".format(nstep_recorded))
		print("{:1.2e} steps after coarsening".format(nstep_coarse))

	### open source trajectory
	with open(datFile, 'r') as fin:
		frames = coarsenFrames(readFrames(fin, nba_total), coarse_time, steps_per_frame, setSteps, 1)

		### check for old coarse file
		nstep_old = 0
		if update and os.path.isfile(outDatFile):
			nstep_old = checkOldFile(outDatFile, frames, nba_total, precision, keepVelocity)
			if not hush:
				print("{:1.2e} steps already coarsened".format(nstep_old))
		else:
			update = False

		### stop before changing anything
		if dryRun:
			if printSteps:
				printFrameSteps(datFile, nba_total, coarse_time, steps_per_frame, setSteps)
			return

		### write coarse trajectory to temporary file
		tmpDatFile = outDatFile + ".tmp"
		with open(tmpDatFile, 'w', buffering=1024*1024) as fout:

			### copy the existing coarse frames
			if update:
				with open(outDatFile, 'r') as fold:
					shutil.copyfileobj(fold, fout)

			### then add the new ones
			for i, (t_line, b_line, e_line, data_lines) in enumerate(frames, nstep_old):

				### keep in touch (status bar only appears if there are new steps)
				if not hush:
					if i == nstep_old:
						initStatusBar("Coarsening trajectory")
						for j in range(nstep_old):
							updateStatusBar(j, nstep_coarse)
					updateStatusBar(i, nstep_coarse)

				### copy header
				fout.write(t_line)
				fout.write(b_line)
				fout.write(e_line)

				### write data
				fout.write(''.join(formatLine(line, precision, keepVelocity) for line in data_lines))

		### swap in the finished file
		os.replace(tmpDatFile, outDatFile)


### check that old coarsened file matches the start of the coarsened trajectory
def checkOldFile(outDatFile, frames, nba_total, precision, keepVelocity):
	nstep_old = 0

	### open old file
	with open(outDatFile, 'r') as fold:

		### loop over old frames
		for t_line_old, _, _, data_old in readFrames(fold, nba_total):
			t_old = int(t_line_old.split()[2])

			### get matching new frame
			frame = next(frames, None)
			if frame is None:
				print(f"Error: Old file has more time steps than new file.\n")
				sys.exit(1)
			t_line_new, _, _, data_new = frame
			t_new = int(t_line_new.split()[2])

			### compare timesteps
			if t_old != t_new:
				print(f"Error: Time step mismatch between old and new file.\n")
				sys.exit(1)

			### get first nucleotide data
			tok_old = data_old[0].split()
			tok_new = formatLine(data_new[0], precision, keepVelocity).split()

			### compare velocity inclusion
			if (len(tok_old) > 9) != (len(tok_new) > 9):
				print(f"Error: Velocity mismatch between old and new file.\n")
				sys.exit(1)

			### compare position
			if tok_old[:3] != tok_new[:3]:
				print(f"Error: Position mismatch between old and new file.\n")
				sys.exit(1)

			### update counter
			nstep_old += 1

	### results
	return nstep_old


### iterate over complete frames
def readFrames(fin, nba_total):
	while True:

		### check for new frame
		t_line = fin.readline()
		if not t_line:
			return

		### read rest of frame
		b_line = fin.readline()
		e_line = fin.readline()
		data_lines = [fin.readline() for _ in range(nba_total)]

		### stop at incomplete frame (simulation still writing)
		if not data_lines[-1].endswith('\n'):
			return

		yield t_line, b_line, e_line, data_lines


### print the step of every frame, next to its step after coarsening
def printFrameSteps(datFile, nba_total, coarse_time, steps_per_frame, setSteps):
	print(f"{'original':<12} coarsened")
	with open(datFile, 'r') as fin:
		for frame, step, step_new in labelFrames(readFrames(fin, nba_total), coarse_time, steps_per_frame, setSteps, 1):
			print(f"{step:<12} {'' if step_new is None else step_new}")


### keep the frames that survive coarsening, optionally renumbering their steps
def coarsenFrames(frames, coarse_time, steps_per_frame, setSteps, offset):
	for (t_line, b_line, e_line, data_lines), step, step_new in labelFrames(frames, coarse_time, steps_per_frame, setSteps, offset):
		if step_new is None:
			continue
		if step_new != step:
			t_line = f"t = {step_new}\n"
		yield t_line, b_line, e_line, data_lines


### pair each frame with its step and its step after coarsening (None if dropped)
def labelFrames(frames, coarse_time, steps_per_frame, setSteps, offset):
	nkept = 0
	for i, frame in enumerate(frames):
		step = int(frame[0].split()[2])

		### mark unused frames
		if coarse_time is not None:
			keep = (i+offset) % coarse_time == 0
		else:
			keep = step % steps_per_frame == 0
		if not keep:
			yield frame, step, None
			continue

		### renumber steps
		step_new = (nkept+offset)*steps_per_frame if setSteps else step
		nkept += 1
		yield frame, step, step_new


### truncate decimals of a nucleotide line, optionally removing velocities
def formatLine(line, precision, keepVelocity):
	parts = line.split() if keepVelocity else line.split(None, 9)[:9]
	out_tokens = []
	for tok in parts:
		dot = tok.find('.')
		if dot == -1:
			out_tokens.append(tok)
			continue
		end = dot + 1 + precision
		e_idx = tok.find('e', dot)
		if e_idx == -1:
			out_tokens.append(tok if end >= len(tok) else tok[:end])
		else:
			out_tokens.append(tok if end >= e_idx else tok[:end] + tok[e_idx:])
	return ' '.join(out_tokens) + '\n'


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

