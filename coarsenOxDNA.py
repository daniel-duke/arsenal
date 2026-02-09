import numpy as np
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

	### command line input
	datFile = sys.argv[1]
	coarse_time = int(sys.argv[2])

	### read old, write new
	readWriteOxDNA(datFile, coarse_time)


################################################################################
### File Managers

### read oxdna trajectory, write coarsened one
def readWriteOxDNA(datFile, coarse_time):
	outDatFile = datFile[:-4] + "_coarse" + datFile[-4:]

	### extract metadata
	print("Getting metadata from trajectory...")
	checkFileExist(datFile, "trajectory", requireData=True)
	with open(datFile, 'r') as f:

		line = f.readline()
		steps_per_record = int(line.split()[2])

		line = f.readline()
		dbox = float(line.split()[2])

		line = f.readline()
		nba_total = 0
		while f.readline().split()[0] != 't':
			nba_total += 1

	### write new file
	print("Writing coarsened trajectory...")
	with open(outDatFile, 'w') as fout:
		with open(datFile) as fin:
			copy = True
			steps_coarse = 0
			steps_recorded = 0
			while True:
				line = fin.readline()
				if not line:
					break
				if len(line.split()) > 1 and line.split()[0] == 't':
					steps_recorded += 1
					if steps_recorded % 100 == 0:
						print("parsed " + str(steps_recorded) + " steps...")
					step = int(line.split()[2])
					if step/steps_per_record % coarse_time == 0:
						copy = True
						steps_coarse += 1
					else:
						copy = False
				if copy == True:
					fout.write(line)

	### report step counts
	print("{:1.2e} steps in simulation".format(steps_recorded*steps_per_record))
	print("{:1.2e} steps recorded".format(steps_recorded))
	print("{:1.2e} steps after coarsening".format(steps_coarse))

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

