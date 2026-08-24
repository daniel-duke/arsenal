import numpy as np
import argparse
import sys
import os

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
	parser.add_argument('coarse_time',type=int)
	parser.add_argument('report',type=int,nargs='?',default=0)
	parser.add_argument('--precision',type=int,default=8)
	parser.add_argument('--keepVelocity',action='store_true')
	args = parser.parse_args()

	### read old, write new
	readWriteOxDNA(args.datFile, args.coarse_time, args.report, args.precision, args.keepVelocity)


################################################################################
### File Managers

### read oxdna trajectory, write coarsened one
def readWriteOxDNA(datFile, coarse_time, report, precision, keepVelocity):
	outDatFile = addSuffix(datFile, "_coarse")

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
		while line[0] != 't':
			nba_total += 1
			line = f.readline()
		steps_per_frame = int(line.split()[2]) - step_initial

	### count lines
	with open(datFile, 'rb') as f:
		nline = sum(1 for _ in f)

	### count steps
	nstep_recorded = nline // (nba_total+3)
	nstep_coarse = nstep_recorded // coarse_time

	### report step counts
	if report: 
		print("{:1.2e} steps in simulation".format(nstep_recorded*steps_per_frame))
		print("{:1.2e} steps in trajectory".format(nstep_recorded))
		print("{:1.2e} steps after coarsening".format(nstep_coarse))
		i = 0

	### write new file
	if report: initStatusBar("Coarsening trajectory")
	with open(datFile, 'r') as fin:
		with open(outDatFile, 'w', buffering=1024*1024) as fout:
			while True:

				### check for new frame
				t_line = fin.readline()
				if not t_line:
					break

				### read rest of frame
				b_line = fin.readline()
				e_line = fin.readline()
				data_lines = [fin.readline() for _ in range(nba_total)]

				### skip unused frames
				step = int(t_line.split()[2])
				if (step // steps_per_frame) % coarse_time != 0:
					continue

				### keep in touch
				if report:
					updateStatusBar(i, nstep_coarse)
					i += 1

				### copy header
				fout.write(t_line)
				fout.write(b_line)
				fout.write(e_line)

				### write data
				chunks = []
				for line in data_lines:
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
					chunks.append(' '.join(out_tokens) + '\n')
				fout.write(''.join(chunks))


################################################################################
### Utility Functions

### determine if file exists
def checkFileExist(file, name="the", required=True, requireData=False):
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
			print(f"Flag: Could not find {name} file.")
			return False


### insert suffix before extension (last dot)
def addSuffix(file, suffix):
	base, ext = os.path.splitext(file)
	return base + suffix + ext


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
		print(f"\n{'='*(length)}")
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
			print("\n", flush=True)


### run the script
if __name__ == "__main__":
	main()
	print()

