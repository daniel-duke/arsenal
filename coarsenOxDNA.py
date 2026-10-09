import numpy as np
import argparse
import sys
import os
import re

## Description
# this script reads an oxDNA trajectory file and rewrites the data into a new
  # file with coarsened time steps.
# this can help reduce file sizes when the output frequency of the simulation
  # was not set appropriately, or when the simulation needed to be run for
  # longer than expected.
# this script loads the old file and writes the new file line by line, which 
  # can be slower for small scripts, but it never times out.
# this script runs from the terminal with two arguments: (1) the trajectory
  # file to coarsen, (2) the coarsening factor.


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('datFile',type=str)
	parser.add_argument('--coarse-time',type=int,default=1)
	parser.add_argument('--precision',type=int,default=8)
	parser.add_argument('--keepVelocity',action='store_true')
	parser.add_argument('--update',action='store_true')
	parser.add_argument('--dryRun',action='store_true')
	parser.add_argument('--hush',action='store_true')
	args = parser.parse_args()

	### read old, write new
	readWriteOxDNA(args.datFile, args.coarse_time, args.hush, args.precision, args.keepVelocity, args.update, args.dryRun)
	if not args.hush: print()


################################################################################
### File Managers

### read oxdna trajectory, write coarsened one
def readWriteOxDNA(datFile, coarse_time, hush, precision, keepVelocity, update, dryRun):
	outDatFile = getCoarseName(datFile, coarse_time)

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
		steps_per_frame = int(line.split()[2]) - step_initial

	### report step counts
	if not hush:
		with open(datFile, 'rb') as f:
			nline = sum(1 for _ in f)
		nstep_recorded = nline // (nba_total+3)
		nstep_coarse = nstep_recorded // coarse_time
		print("{:1.2e} steps in simulation".format(nstep_recorded*steps_per_frame))
		print("{:1.2e} steps in trajectory".format(nstep_recorded))
		print("{:1.2e} steps after coarsening".format(nstep_coarse))

	### open source trajectory
	with open(datFile, 'r') as fin:
		frames = readFrames(fin, nba_total, steps_per_frame, coarse_time)

		### check for old coarse file
		nstep_old = 0
		if update and os.path.isfile(outDatFile):
			nstep_old, nbyte_old = checkOldFile(outDatFile, frames, nba_total, precision, keepVelocity)
			if not hush:
				print("{:1.2e} steps already coarsened".format(nstep_old))
		else:
			update = False

		### stop before changing anything
		if dryRun:
			return

		### remove any incomplete frame from old file
		if update:
			os.truncate(outDatFile, nbyte_old)

		### write coarse trajectory
		with open(outDatFile, 'a' if update else 'w', buffering=1024*1024) as fout:
			for i, (t_line, b_line, e_line, data_lines) in enumerate(frames, nstep_old):

				### keep in touch (status bar only appears if there are new steps)
				if not hush:
					if i == nstep_old:
						initStatusBar("Coarsening trajectory")
					updateStatusBar(i, nstep_coarse)

				### copy header
				fout.write(t_line)
				fout.write(b_line)
				fout.write(e_line)

				### write data
				fout.write(''.join(formatLine(line, precision, keepVelocity) for line in data_lines))


### check that old coarsened file matches the start of the coarsened trajectory
def checkOldFile(outDatFile, frames, nba_total, precision, keepVelocity):
	nstep_old = 0
	nbyte_old = 0

	### open old file
	with open(outDatFile, 'r') as fold:

		### loop over old frames
		for t_line_old, _, _, data_old in readFrames(fold, nba_total, 1, 1):
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

			### update counters
			nstep_old += 1
			nbyte_old = fold.tell()

	### results
	return nstep_old, nbyte_old


### iterate over complete frames that survive coarsening
def readFrames(fin, nba_total, steps_per_frame, coarse_time):
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

		### skip unused frames
		step = int(t_line.split()[2])
		if (step // steps_per_frame) % coarse_time != 0:
			continue

		yield t_line, b_line, e_line, data_lines


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
			sys.exit()
	else:
		if required:
			print(f"Error: Could not find {name} file:")
			print(file + "\n")
			sys.exit()
		else:
			if not hush:
				print(f"Flag: Could not find {name} file.")
			return False


### name coarse file, tracking coarsening factors in the suffix
def getCoarseName(file, coarse_time):
	base, ext = os.path.splitext(file)
	_, sep, suffix = os.path.basename(base).rpartition('_')
	if sep and re.fullmatch(r'cx\d+(x\d+)*', suffix):
		return f"{base}x{coarse_time}{ext}"
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

