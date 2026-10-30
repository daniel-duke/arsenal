import arsenal as ars
import numpy as np
import argparse
import sys

## Description
# this script reads an oxDNA topology (.top) and configuration (.dat) file and
  # writes a LAMMPS geometry file (.in) for simulating the structure in LAMMPS.
# the geometry is written in oxDNA units, with one ellipsoid per nucleotide.
# atom types are set by base (A=1, C=2, G=3, T=4), molecule indices are the
  # strand numbers, and each nucleotide is bonded to its 5' neighbor.
# the output is identical to what tacoxDNA.py writes with the --lammps option.
# note that this differs from lammpify.py, which is meant for visualization.
# this script runs from the command line with two arguments: (1) the topology
  # file, (2) the configuration file.
# the script outputs a file with the same name as the configuration file, but
  # with the extension ".in".


################################################################################
### Parameters

### lammps oxDNA constants (oxDNA units)
mass_nt = 3.1575
diam_nt = 1.173984503142341
base_types = {'A':1, 'C':2, 'G':3, 'T':4}


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('topFile',	type=str,	help="name of topology file")
	parser.add_argument('datFile',	type=str,	help="name of configuration file")
	args = parser.parse_args()

	### read topology and configuration
	strands, bases, n5s = readTop(args.topFile)
	points, dbox, axes = ars.readOxDNA(args.datFile, getAxes=True, keepOxUnits=True, ignorePBC=True)
	r = points[0]
	axes = axes[0]

	### check consistency
	nba_total = len(strands)
	if len(r) != nba_total:
		print("Error: Topology and configuration files contain different numbers of bases.\n")
		sys.exit(1)

	### set per-nucleotide properties
	types = np.array([ base_types[b] for b in bases ])
	quats = np.array([ ars.axesToQuat(a) for a in axes ])
	quats = np.round(quats, 8) + 0.0
	radii = np.full((nba_total,3), diam_nt/2)

	### set bonds (to 5' neighbor, so each bond is only listed once)
	bonds = np.array([ [1,i+1,n5s[i]+1] for i in range(nba_total) if n5s[i] >= 0 ])

	### write geometry
	print("Writing geometry...")
	geoFile = ars.changeExtension(args.datFile, "in")
	ars.writeGeo(geoFile, dbox, r, strands, types, bonds, natomType=4, masses=[mass_nt]*4, radii=radii, quats=quats, x_precision=8)


################################################################################
### File Handlers

### read oxDNA topology file
def readTop(topFile):

	### notes
	# strands are 1-based, as they appear in the topology file.
	# n5s gives the (0-based) index of each base's 5' neighbor, or -1 if there is none.

	### read file
	ars.checkFileExist(topFile, "topology")
	with open(topFile, 'r') as f:
		content = f.readlines()

	### parse header
	nba_total = int(content[0].split()[0])

	### parse bases
	strands = np.zeros(nba_total, dtype=int)
	bases = [None]*nba_total
	n5s = np.zeros(nba_total, dtype=int)
	for i in range(nba_total):
		line = content[i+1].split()
		strands[i] = int(line[0])
		bases[i] = line[1]
		n5s[i] = int(line[3])

	### result
	return strands, bases, n5s


### run the script
if __name__ == "__main__":
	main()
	print()
