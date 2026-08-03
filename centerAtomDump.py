import numpy as np
import argparse
import sys
import os

## Description
# this script reads a lammps-style atom dump, centers the atoms (either
  # by atom com, molecule com, or given molecule ID), optionally unwraps
  # the atoms from the boundary by molecule, and optionally sets the color
  # in the new trajectory.
# the script outputs a file with the same name as the input file, modified to
  # with "_centered" at the end (but before ".dat"); for example, the input file
  # "trajectory.dat" would result in output file "trajectory_centered.dat"


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('datFile',type=str)
	parser.add_argument('center',type=str,nargs='?',default='com')
	parser.add_argument('unwrap',type=int,nargs='?',default=0)
	parser.add_argument('setColor',type=int,nargs='?',default=0)
	parser.add_argument('report',type=int,nargs='?',default=0)
	parser.add_argument('--mol-col',type=int,default=2)
	parser.add_argument('--color-col',type=int,default=2)
	parser.add_argument('--excludeDummy',action='store_true')
	args = parser.parse_args()

	### convert to 0-based column indices
	col_mol = args.mol_col-1
	col_color = args.color_col-1

	### read file (no skipped steps, no coarsening)
	points, molecules, colors, dbox3s, steps_per_frame = readAtomDump(args.datFile, col_mol, col_color, args.report)

	### center the points
	points_centered = centerPointsMolecule(points, molecules, dbox3s, args.center, args.unwrap, args.report, args.excludeDummy)

	### write output file
	outDatFile = addSuffix(args.datFile, "_centered")
	writeAtomDump(outDatFile, points_centered, colors, dbox3s, steps_per_frame, args.setColor, args.report)


################################################################################
### File Managers

### read lammps-style atom dump
def readAtomDump(datFile, col_mol, col_color, report):

	### count lines
	checkFileExist(datFile, "trajectory")
	with open(datFile, 'rb') as f:
		nline = sum(1 for _ in f)

	### extract metadata from header
	with open(datFile, 'r') as f:

		### read header
		header = []
		for i in range(9):
			header.append(f.readline())
		step_init = int(header[1].split()[0])
		nbd_total = int(header[3].split()[0])

		### parse columns
		col_x = None
		line = header[8].split()
		for i in range(2,len(line)):
			if line[i] == 'xs':
				isScaled = True
				col_x = i-2
			elif line[i] == 'x':
				isScaled = False
				col_x = i-2
		if col_x is None:
			print("Error: No position data found.")
			sys.exit()

		### look for next frame
		steps_per_frame = 0
		for i in range(nbd_total+2):
			line = f.readline()
		if line: steps_per_frame = int(line.split()[0]) - step_init

	### count steps
	nstep = nline // (nbd_total+9)

	### initialize
	points = np.zeros((nstep,nbd_total,3))
	molecules = np.zeros(nbd_total,dtype=int)
	colors = np.zeros(nbd_total,dtype=int)
	dbox3s = np.zeros((nstep,3))

	### extract the data
	if report: initStatusBar("Reading file")
	with open(datFile, 'r') as f:
		for i in range(nstep):
			content = [f.readline() for _ in range(nbd_total+9)]
			for k in range(3):
				line = content[5+k].split()
				dbox3s[i,k] = float(line[1]) - float(line[0])
			for j in range(nbd_total):
				line = content[9+j].split()
				if i == 0:
					molecules[j] = int(line[col_mol])
					colors[j] = int(line[col_color])
				points[i,j] = np.array(line[col_x:col_x+3],dtype=float)
				if isScaled:
					points[i,j] = (points[i,j]-1/2)*dbox3s[i]
			points[i] = applyPBC(points[i], dbox3s[i])
			if report: updateStatusBar(i,nstep)

	### result
	return points, molecules, colors, dbox3s, steps_per_frame


### write lammps-style trajectory
def writeAtomDump(outDatFile, points, col2s, dbox3s, steps_per_frame, setColor, report):
	nstep = points.shape[0]
	npoint = points.shape[1]
	len_npoint = len(str(npoint))
	len_ncol2 = len(str(max(col2s)))
	if report: initStatusBar("Writing file")
	with open(outDatFile,'w') as f:
		for i in range(nstep):
			len_dbox = len(str(int(max(dbox3s[i])/2)))
			f.write(f"ITEM: TIMESTEP\n{i*steps_per_frame}\n")
			f.write(f"ITEM: NUMBER OF ATOMS\n{npoint}\n")
			f.write(f"ITEM: BOX BOUNDS pp pp pp\n")
			f.write(f"-{dbox3s[i,0]/2:<{len_dbox+3}.2f} {dbox3s[i,0]/2:<{len_dbox+3}.2f} xlo xhi\n")
			f.write(f"-{dbox3s[i,1]/2:<{len_dbox+3}.2f} {dbox3s[i,1]/2:<{len_dbox+3}.2f} ylo yhi\n")
			f.write(f"-{dbox3s[i,2]/2:<{len_dbox+3}.2f} {dbox3s[i,2]/2:<{len_dbox+3}.2f} zlo zhi\n")
			if setColor:
				f.write("ITEM: ATOMS id type xs ys zs\n")
			else:
				f.write("ITEM: ATOMS id mol xs ys zs\n")
			for j in range(npoint):
				f.write(f"{j+1:<{len_npoint}} " + \
						f"{col2s[j]:<{len_ncol2}}  " + \
						f"{points[i,j,0]/dbox3s[i,0]+1/2:10.8f} " + \
						f"{points[i,j,1]/dbox3s[i,1]+1/2:10.8f} " + \
						f"{points[i,j,2]/dbox3s[i,2]+1/2:10.8f}\n")
			if report: updateStatusBar(i,nstep)


################################################################################
### Calculations

### shift trajectory, placing the given point at the center, optionally unwrapping molecules at boundary
def centerPointsMolecule(points, molecules, dbox3s, center, unwrap, report, excludeDummy):
	nstep = points.shape[0]
	npoint = points.shape[1]
	nmolecule = int(max(molecules))

	### sort points by molecule
	points_moleculed = sortPointsByMolecule(points, molecules)

	### initialize
	molecule_coms = np.zeros((nmolecule,3))
	points_centered = np.zeros((nstep,npoint,3))

	### loop over steps
	if report: initStatusBar("Centering")
	for i in range(nstep):

		### calculate molecule coms
		for j in range(nmolecule):
			if checkAnyDummy(points_moleculed[j][i]):
				if not excludeDummy:
					print("Warning: Potential dummy beads detected.")
				elif not checkAllDummy(points_moleculed[j][i]):
					print("Error: Molecule contains mixed dummy and activated beads.\n")
					sys.exit()
			molecule_coms[j] = calcCOM(points_moleculed[j][i], dbox3s[i], excludeDummy)

		### set centering point
		if center == 'none':
			com = np.zeros(3)
		elif center == 'com' or center == 'com_points' or center == 'com_beads' or center == 'com_bases':
			com = calcCOM(points[i], dbox3s[i], excludeDummy)
		elif center == 'com_molecules' or center == 'com_clusters':
			com = calcCOM(molecule_coms, dbox3s[i], excludeDummy)
		elif isinteger(center) and int(center) <= nmolecule:
			com = molecule_coms[int(center)-1]
		else:
			print("Error: Cannot center points - center must be either 'none', 'com', 'com_molecules', or integer <= nmolecule.\n")
			sys.exit()

		### center the points
		per_point_coms = molecule_coms[molecules-1]
		dummy = np.all(per_point_coms == 0, axis=1)
		points_centered[i,~dummy] = applyPBC(points[i,~dummy] - com, dbox3s[i])

		### unwrap molecules at boundary
		if unwrap:
			molecule_coms_centered = np.zeros((nmolecule,3))
			for j in range(nmolecule):
				if not all(molecule_coms[j]==0):
					molecule_coms_centered[j] = applyPBC(molecule_coms[j]-com, dbox3s[i])
			refs = molecule_coms_centered[molecules-1]
			points_centered[i] = refs + applyPBC(points_centered[i]-refs, dbox3s[i])

		### progress updates
		if report: updateStatusBar(i,nstep)

	### result
	return points_centered


### calculate center of mass, using method from Bai and Breen 2008
def calcCOM(r, dbox3, excludeDummy=False):
	if excludeDummy and len(r.shape)>1: r = r[~np.all(r==0, axis=1)]
	xi_bar = np.mean( np.cos(2*np.pi*(r/dbox3+1/2)), axis=0 )
	zeta_bar = np.mean( np.sin(2*np.pi*(r/dbox3+1/2)), axis=0 )
	theta_bar = np.arctan2(-zeta_bar, -xi_bar) + np.pi
	r_ref = dbox3*(theta_bar/(2*np.pi)-1/2)
	com = r_ref + np.mean( applyPBC(r-r_ref, dbox3), axis=0 )
	return com


################################################################################
### Utility Functions

### apply periodic boundary condition
def applyPBC(r, dbox):
	return r - dbox*np.round(r/dbox)


### sort points into molecules
def sortPointsByMolecule(points, molecules):
	nstep = points.shape[0]
	npoint = points.shape[1]
	nmolecule = int(max(molecules))
	points_moleculed = [None]*nmolecule
	for m in range(nmolecule):
		points_moleculed[m] = np.zeros((nstep,sum(molecules==m+1),3))
	for i in range(nstep):
		n_molecule_count = np.zeros(nmolecule, dtype=int)
		for j in range(npoint):
			points_moleculed[molecules[j]-1][i,n_molecule_count[molecules[j]-1]] = points[i,j]
			n_molecule_count[molecules[j]-1] += 1
	return points_moleculed


### determine if position array contains any dummy beads
def checkAnyDummy(r):
	return np.any(np.all(r==0, axis=1))


### determine if position array contains all dummy beads
def checkAllDummy(r):
	return np.all(r==0)


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


### test if variable is an integer (python int, numpy int, and string int all count)
def isinteger(x):
	if isinstance(x, int) or isinstance(x, np.int64):
		return True
	if isinstance(x, str):
		try:
			int(x)
			return True
		except ValueError:
			return False
	return False


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

