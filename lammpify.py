import arsenal as ars
import numpy as np
from scipy import stats
import argparse
import shutil
import re
import sys
import os

## Description
# this script reads an oxDNA or lammps trajectory file, discards unnecessary
  # information (for oxDNA), fixes any boundary jumps, centers the positions,
  # recolors the particles, and rewites the data in LAMMPS format.
# cluster files contain lists of base indices, with one list every other line,
  # generally resembling something like "CLUSTER 1 \n 1 2 3 \n CLUSTER 2..."
# to create a cluster file for oxDNA strucutres, select the bases in oxView and
  # select "Download Selection IDs" then add the list of bases to your cluster
  # file, which should be formatted as described above.
# to create a cluster file for LAMMPS structures, load the file in OVITO and 
  # add the "Particle Identifier" for the desired particles to the cluster file.
# note that oxDNA base indices are 0 based and LAMMPS bead indices are 1 based.


################################################################################
### Parameters

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--geoFile',		type=str, default=None,		help="for lammps simulations, name of geometry file")
	parser.add_argument('--datFile',		type=str, default=None,		help='for lammps simulations, name of trajectory file')
	parser.add_argument('--oxFiles',		type=str, nargs=2,			help='for oxdna simulations, name of topology and configuration/trajectory files')
	parser.add_argument('--clusterFile',	type=str, default=None,		help="name of clusters file")
	parser.add_argument('--forceFile',		type=str, default=None,		help="name of external forces file")
	parser.add_argument('--nstep-skip',		type=int, default=0, 		help="number of recorded initial steps to skip")
	parser.add_argument('--coarse-time',	type=int, default=1, 		help="coarse factor for time steps")
	parser.add_argument('--center',			action='store_true',		help="whether to center the trajectory")
	parser.add_argument('--align',			action='store_true',		help="whether to align principal components with coordinate axes (first frame only)")
	parser.add_argument('--split',			action='store_true',		help="whether split oxDNA particles into backbone and base beads")
	parser.add_argument('--style',			type=str, default='mol', 	help="for lammps simulations, atom style in geometry file (mol, full, ox)")
	parser.add_argument('--units',			type=str, default='nm', 	help="for lammps simuations, length scale of input data (nm, ang, ox)")
	parser.add_argument('--coloring',		type=str, default='scaf', 	help="for oxDNA simulations, coloring (scaf, strand, seq)")
	parser.add_argument('--dbox',			type=int, default=None, 	help="new box diameter (not used for centering)")

	### set arguments
	args = parser.parse_args()
	geoFile = args.geoFile
	datFile = args.datFile
	oxFiles = args.oxFiles
	clusterFile = args.clusterFile
	forceFile = args.forceFile
	nstep_skip = args.nstep_skip
	coarse_time = args.coarse_time
	center = args.center
	align = args.align
	style = args.style
	units = args.units
	split = args.split
	coloring = args.coloring
	dbox_new = args.dbox

	### check for conflicting inputs
	if geoFile is not None and oxFiles is not None:
		print("Error: Must provide either geometry file or oxDNA files, not both.")
		sys.exit()

	### adjustments for oxDNA simulations
	if oxFiles is not None:
		split = True

	### adjustments for split lammps simulation
	if geoFile is not None and split: 
		style = 'ox'
		units = 'ox'


################################################################################
### Main

	### output file parameters
	outFold = "visualize/"
	outGeoFile = outFold + "geometry.in"
	outDatFile = outFold + "trajectory.dat"

	### create output folder
	ars.createEmptyFold(outFold)

	### lammpify lammps
	if geoFile is not None:

		### read geometry
		if style == 'ox':
			points_init, dbox3, molecules, types, radii, quats_init, bonds = ars.readGeo(geoFile, style=style, getDbox3=True, getBonds=True)
		else:
			points_init, dbox3, molecules, types, bonds = ars.readGeo(geoFile, style=style, getDbox3=True, getBonds=True)

		### read trajectory
		if datFile is not None:
			if not split:
				points, dbox3s = ars.readAtomDump(datFile, nstep_skip, coarse_time, getDbox3s=True)
			else:
				points, dbox3s, quats = ars.readAtomDump(datFile, nstep_skip, coarse_time, getDbox3s=True, getQuats=True)

		### unit conversion
		scale = convertToNm(units)
		points_init *= scale
		dbox3 *= scale
		if datFile is not None:
			points *= scale
			dbox3s *= scale

		### center trajectory
		if center:
			points_init = ars.centerPointsMolecule(points_init, molecules, dbox3, center='com', unwrap=False)
			if datFile is not None:
				points = ars.centerPointsMolecule(points, molecules, dbox3s, center='com', unwrap=False)

		### align first frame
		if align:
			points_init, PCs_init = ars.alignPCs(points_init, getPCs=True)
			if datFile is not None:
				points, PCs = ars.alignPCs(points, getPCs=True)

		### set colors
		if clusterFile is not None:
			clusters = ars.readCluster(clusterFile)
			colors = getMoleculesFromClustersB1(clusters, len(types))
		else:
			colors = types

		### set box diameter
		if dbox_new is not None:
			dbox3 = dbox_new

		### beads
		if not split:

			### write output
			ars.writeGeo(outGeoFile, dbox3, points_init, molecules, colors, bonds)
			if datFile is not None:
				ars.writeAtomDump(outDatFile, dbox3s, points, colors)

		### nucleotides
		else:

			### split and write geometry
			axes_init = quatsToAxes(quats_init)
			if align: axes_init = axes_init @ PCs
			points_init, quats_init, molecules, colors, radii, bonds = splitNucleotides(points_init, axes_init, molecules, colors, bonds)
			ars.writeGeo(outGeoFile, dbox3, points_init, molecules, colors, bonds, radii=radii, quats=quats_init)

			### split and write trajectory
			if datFile is not None:
				axes = quatsToAxes(quats)
				points, quats = splitNucleotides(points, axes, molecules, colors, bonds)[:2]
				ars.writeAtomDump(outDatFile, dbox3s, points, colors, quats=quats)

	### lammpify oxdna
	elif oxFiles is not None:

		### read topology
		strands, bases, bonds, nba_total = readTop(oxFiles[0])

		### read trajecotry
		points, dbox3, axes = ars.readOxDNA(oxFiles[1], nstep_skip, coarse_time, getDbox3=True, getAxes=True)
		nstep = points.shape[0]

		### center trajectory
		if center: points = ars.centerPointsMolecule(points, strands, dbox3, center='com', unwrap=True)

		### align first frame
		if align:
			points, PCs = ars.alignPCs(points, getPCs=True)
			if split:
				axes = axes @ PCs[:,None]

		### set colors
		if clusterFile is not None:
			clusters = ars.readCluster(clusterFile)
			colors = getMoleculesFromClustersB0(clusters, nba_total)
		elif forceFile is not None:
			clusters = readForce(forceFile)
			colors = getMoleculesFromClustersB0(clusters, nba_total)
		elif coloring == 'scaf':
			strand_scaffold = stats.mode(strands).mode
			colors = np.where(strands == strand_scaffold,1,2)
		elif coloring == 'strand':
			colors = strands
		elif coloring == 'seq':
			colors = bases
		else:
			print("Error: Unrecognized coloring.")
			sys.exit()

		### set box diameter
		if dbox_new is not None:
			dbox3 = dbox_new

		### beads
		if not split:

			### write output
			ars.writeGeo(outGeoFile, dbox3, points[0], strands, colors, bonds)
			if points.shape[0] > 1:
				ars.writeAtomDump(outDatFile, dbox3, points, colors)

		### nucleotides
		else:

			### split data
			points, quats, strands, colors, radii, bonds = splitNucleotides(points, axes, strands, colors, bonds)

			### write output
			ars.writeGeo(outGeoFile, dbox3, points[0], strands, colors, bonds, radii=radii, quats=quats[0])
			if nstep > 1:
				ars.writeAtomDump(outDatFile, dbox3, points, colors, quats=quats)

	### error
	else:
		print("Error: Unknown simulation type, try again.")
		sys.exit()


################################################################################
### File Handlers

### read oxdna topology
def readTop(topFile):
	base_types = {'A':1, 'C':2, 'G':3, 'T':4}
	ars.checkFileExist(topFile, "topology")
	with open(topFile) as f:
		content = f.readlines()
	nba_total = int(content[0].split()[0])
	strands = np.zeros(nba_total,dtype=int)
	bases = np.zeros(nba_total,dtype=int)
	bonds = np.ones((nba_total,3),dtype=int)
	bond_count = 0
	for i in range(nba_total):
		line = content[i+1].split()
		strands[i] = int(line[0])
		bases[i] = base_types[line[1]]
		if int(content[i+1].split()[3]) != -1:
			bonds[bond_count,1] = i+1
			bonds[bond_count,2] = int(line[3])+1
			bond_count += 1
	bonds = bonds[:bond_count,:]
	return strands, bases, bonds, nba_total


################################################################################
### Utility Functions

### unit conversion factor to nm
def convertToNm(units):

	### parse unit type
	if units == 'nm':
		scale = 1
	elif units == 'ox':
		scale = ars.ox2nm
	elif units == 'ang':
		scale = 0.1
	else:
		print("Error: Unrecognized input data units.")
		sys.exit()

	### result
	return scale


### read COM force file to get clusters
def readForce(forceFile):
	ars.checkFileExist(forceFile, "force")
	with open(forceFile, 'r') as f:
		content = f.read()
	
	### split content
	forces = re.findall(r'\{(.*?)\}', content, re.DOTALL)

	### loop over forces
	clusters = []
	for force in forces:
		com_line = re.search(r'com_list\s*=\s*([0-9,\s]+)', force)
		ref_line = re.search(r'ref_list\s*=\s*([0-9,\s]+)', force)

		if com_line is None or ref_line is None:
			continue

		com_bais = [int(i) for i in com_line.group(1).replace(',', ' ').split()]
		ref_bais = [int(i) for i in ref_line.group(1).replace(',', ' ').split()]

		clusters.append(com_bais)
		clusters.append(ref_bais)

	### results
	return clusters


### split particles into backbone and base sites, adding appropriate bonds
def splitNucleotides(points, axes, molecules, colors, bonds):

	### add time dimension to single frames
	points, ndim_add = ars.padDims(points,3)
	axes = ars.padDims(axes,4)[0]

	### counts
	nstep = points.shape[0]
	npoint = points.shape[1]

	### initialize
	points_split = np.zeros((nstep,npoint*2,3))
	quats = np.zeros((nstep,npoint*2,4))
	quats[:,:,0] = 1

	### split position data
	points_split[:,:npoint] = points + ars.ox2nm*(-0.34*axes[:,:,0] + 0.3408*axes[:,:,1])
	points_split[:,npoint:] = points + ars.ox2nm*(0.4*axes[:,:,0])

	### calculate orientation data
	ars.initStatusBar("Calculating quaternions")
	for i in range(nstep):
		for j in range(npoint):
			quats[i,npoint+j] = ars.axesToQuat(axes[i,j])
		ars.updateStatusBar(i,nstep)

	### add base to backbone bonds
	for j in range(npoint):
		bonds = np.append(bonds,[[2,j+1,npoint+j+1]], axis=0)

	### identification data
	molecules = np.concatenate((molecules,molecules))
	colors = np.concatenate((colors,colors))

	### bead sizes
	radii = np.full((npoint*2,3), 0.34)
	radii[npoint:,2] = 0.17

	### add time dimension to single frames
	points_split = ars.trimDims(points_split, ndim_add)
	quats = ars.trimDims(quats, ndim_add)

	### results
	return points_split, quats, molecules, colors, radii, bonds


### wrapper function for converting quaternion trajectory to coordinate axes trajectory
def quatsToAxes(quats):

	### add time dimension to single frames
	quats, ndim_add = ars.padDims(quats)

	### initialize
	nstep = quats.shape[0]
	npoint = quats.shape[1]
	axes = np.zeros((nstep,npoint,3,3))

	### calculations
	ars.initStatusBar("Calculating axes")
	for i in range(nstep):
		for j in range(npoint):
			axes[i,j] = ars.quatToAxes(quats[i,j])
		ars.updateStatusBar(i,nstep)

	### add time dimension to single frames
	quats = ars.trimDims(quats, ndim_add)

	### result
	return axes


### identify molecules (1 for unidentified, 2+ for molecule IDs)
def getMoleculesFromClustersB0(clusters, npoint):
	molecules = np.ones(npoint,dtype=int)
	for c in range(len(clusters)):
		for j in range(len(clusters[c])):
			index = clusters[c][j]
			if index >= npoint:
				print("Error: requested index " + str(index) + " exceeds the number of beads in the simulation (" + str(npoint) + ").")
				sys.exit()
			molecules[index] = c+2
	return molecules


### identify molecules (1 for unidentified, 2+ for molecule IDs)
def getMoleculesFromClustersB1(clusters, npoint):
	molecules = np.ones(npoint,dtype=int)
	for c in range(len(clusters)):
		for j in range(len(clusters[c])):
			index = clusters[c][j]
			if index > npoint:
				print("Error: requested index " + str(index) + " exceeds the number of beads in the simulation (" + str(npoint) + ").")
				sys.exit()
			molecules[index-1] = c+2
	return molecules


### run the script
if __name__ == "__main__":
	main()
	print()

