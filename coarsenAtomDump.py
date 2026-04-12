import numpy as np
import sys
import os

## Description
# this script reads a LAMMPS-style trajectory file and rewrites the data into a
  # new file with coarsened time steps, as set by the coarse_time parameter; 
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

	### command line input
	datFile = sys.argv[1]
	coarse_time = int(sys.argv[2])

	### optional third argument
	report_every = 0
	if len(sys.argv) > 3:
		report_every = int(sys.argv[3])

	### read old, write new
	readWriteAtomDump(datFile, coarse_time, report_every)


################################################################################
### File Managers

### read LAMMPS-style atom dump, write coarsened one
def readWriteAtomDump(datFile, coarse_time, report_every):
	outDatFile = datFile[:-4] + "_coarse" + datFile[-4:]

	### extract metadata
	if report_every: print("Getting metadata from trajectory...")
	checkFileExist(datFile, "trajectory", requireData=True)
	with open(datFile, 'r') as f:

		for i in range(2): line = f.readline()
		step_initial = int(line.split()[0])

		for i in range(2): line = f.readline()
		nbd_total = int(line.split()[0])

		for i in range(2): line = f.readline()
		dbox = float(line.split()[1])

		for i in range(nbd_total+5): line = f.readline()
		steps_per_frame = int(line.split()[0]) - step_initial

	### write new file
	if report_every: print("Writing coarsened trajectory...")
	with open(outDatFile, 'w') as fout:
		with open(datFile) as fin:
			copy = True
			step_count = 0
			step_count_coarse = 0
			while True:
				line = fin.readline()
				if not line:
					break
				if len(line.split()) > 1 and line.split()[1] == 'TIMESTEP':
					step_count += 1
					if report_every and step_count % report_every == 0:
						print("parsed " + str(step_count) + " steps...")
					current_position = fin.tell()
					next_line = fin.readline()
					step = int(next_line.split()[0])
					fin.seek(current_position)
					if step/steps_per_frame % coarse_time == 0:
						copy = True
						step_count_coarse += 1
					else:
						copy = False
				if copy == True:
					fout.write(line)

	### report step counts
	if report_every: 
		print("{:1.2e} steps in simulation".format(step_count*steps_per_frame))
		print("{:1.2e} steps in trajectory".format(step_count))
		print("{:1.2e} steps after coarsening".format(step_count_coarse))

	### result
	return


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


### run the script
if __name__ == "__main__":
	main()
	print()

