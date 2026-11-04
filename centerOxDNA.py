import numpy as np
import argparse
import sys
import os

## Description
# this script reads an oxDNA trajectory, centers the bases, and writes the data
  # into a new trajectory file.
# this script runs from the terminal with the trajectory file to center as its
  # argument; the number of decimal places written is set with --precision
  # (at most 15, the same as what oxDNA prints).
# the output goes next to the input (with "_centered" added to the name)
  # unless --out is given.
# the output is written to a temporary file and then swapped in, so anything
  # reading the output (e.g. rsync) never sees a half-written file.


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('datFile',type=str)
	parser.add_argument('--out',type=str,default=None)
	parser.add_argument('--precision',type=int,default=8)
	parser.add_argument('--hush',action='store_true')
	args = parser.parse_args()

	### check arguments
	if args.precision > 15:
		print("Error: Precision should not be greater than 15.\n")
		sys.exit(1)

	### read old, write new
	readWriteOxDNA(args.datFile, args.out, args.precision, args.hush)
	if not args.hush: print()


################################################################################
### File Managers

### read oxdna trajectory, write centered one
def readWriteOxDNA(datFile, outDatFile, precision, hush):
	if outDatFile is None:
		outDatFile = addSuffix(datFile, "_centered")

	### extract metadata
	checkFileExist(datFile, "trajectory", requireData=True)
	with open(datFile, 'r') as f:

		for i in range(2): line = f.readline()
		dbox3 = np.array(line.split()[2:5],dtype=float)

		line = f.readline()
		line = f.readline()
		ncol = len(line.split())
		nba_total = 0
		while line and line[0] != 't':
			nba_total += 1
			line = f.readline()

	### count lines
	with open(datFile, 'rb') as f:
		nline = sum(1 for _ in f)

	### count steps
	nstep = nline // (nba_total+3)

	### write new file to temporary file
	data = np.zeros((nba_total,ncol),dtype=float)
	if not hush: initStatusBar("Centering")
	tmpDatFile = outDatFile + ".tmp"
	with open(tmpDatFile, 'w', buffering=1024*1024) as fout:
		with open(datFile) as fin:
			for i in range(nstep):
				for j in range(3):
					line = fin.readline()
					if not line:
						break
					fout.write(line)
				if not line:
					break
				for j in range(nba_total):
					data[j] = fin.readline().split()
				com = calcCOM(data[:,:3], dbox3)
				data[:,:3] = applyPBC( data[:,:3]-com, dbox3 )
				np.savetxt(fout, data, fmt=f'%.{precision}f')
				if not hush: updateStatusBar(i,nstep)

	### swap in the finished file
	os.replace(tmpDatFile, outDatFile)


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

