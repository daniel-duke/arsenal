import numpy as np
import sys
import os

## Description
# this script reads an oxDNA trajecotry, centers the bases, and writes the data
  # into a new trajectory file.
# this script runs from the command line with one argumnets, the name of
  # the trajectory file to load and center.
# there is also an optional second argument which sets the precision to use when
  # writing the output file (defaults to 15, the same as what oxDNA prints).
# there is also an optional third argument which sets the number of steps
  # between printing progress updates (defaults to 0 for no updates).
# the script outputs a file with the same name as the input file, modified to
  # with "_centered" at the end (but before ".dat"); for example, the input file
  # "trajectory.dat" would result in output file "trajectory_centered.dat"


################################################################################
### Heart

def main():

	### command line input
	datFile = sys.argv[1]

	### optional second argument
	precision = 15
	if len(sys.argv) > 2:
		precision = int(sys.argv[2])
		if precision > 15:
			print("Error: Precision should not be greater than 15.")
			sys.exit()

	### optional third argument
	report_every = 0
	if len(sys.argv) > 3:
		report_every = int(sys.argv[3])

	### read old, write new
	readWriteOxDNA(datFile, precision, report_every)


################################################################################
### File Managers

### read oxdna trajectory
def readWriteOxDNA(datFile, precision, report_every):
	outDatFile = datFile[:-4] + "_centered" + datFile[-4:]

	### extract metadata
	if report_every: print("Getting metadata from trajectory...")
	checkFileExist(datFile, "trajectory")
	with open(datFile, 'r') as f:

		for i in range(2): line = f.readline()
		dbox3 = np.array(line.split()[2:5],dtype=float)

		line = f.readline()
		nba_total = 0
		while f.readline().split()[0] != 't':
			nba_total += 1

	### write new file
	if report_every: print("Writing centered trajectory...")
	data = np.zeros((nba_total,15),dtype=float)
	step_count = 0
	with open(outDatFile, 'w') as fout:
		with open(datFile) as fin:
			while True:
				for i in range(3):
					line = fin.readline()
					if not line:
						return
					fout.write(line)
				step_count += 1
				for j in range(nba_total):
					data[j] = np.fromstring(fin.readline(), sep=' ')
				com = calcCOM(data[:,:3], dbox3)
				data[:,:3] = applyPBC( data[:,:3]-com, dbox3 )
				np.savetxt(fout, data, fmt=f'%.{precision}g')
				if report_every and step_count % report_every == 0:
					print("centered " + str(step_count) + " steps...")


################################################################################
### Utility Functions

### calculate center of mass, using method from Bai and Breen 2008
def calcCOM(r, dbox3):
	xi_bar = np.mean( np.cos(2*np.pi*(r/dbox3+1/2)), axis=0 )
	zeta_bar = np.mean( np.sin(2*np.pi*(r/dbox3+1/2)), axis=0 )
	theta_bar = np.arctan2(-zeta_bar, -xi_bar) + np.pi
	r_ref = dbox3*(theta_bar/(2*np.pi)-1/2)
	com = r_ref + np.mean( applyPBC(r-r_ref, dbox3), axis=0 )
	return com


### apply periodic boundary condition
def applyPBC(r, dbox):
	return r - dbox*np.round(r/dbox)


### test if files exist
def checkFileExist(file, name="the", required=True):
	if os.path.isfile(file):
		return True
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


