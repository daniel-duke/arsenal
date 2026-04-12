import arsenal as ars
import numpy as np
from scipy import stats
import argparse
import shutil
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
	parser.add_argument('--oxFiles',		type=str, nargs=2,			help='for oxdna simulations, name of topology and configuration files')
	parser.add_argument('--clusterFile',	type=str, default=None,		help="name of clusters file")
	parser.add_argument('--nstep-skip',		type=int, default=0, 		help="number of recorded initial steps to skip")
	parser.add_argument('--coarse-time',	type=int, default=1, 		help="coarse factor for time steps")
	parser.add_argument('--center',			action='store_true',		help="whether to center the trajectory")

	### set arguments
	args = parser.parse_args()
	geoFile = args.geoFile
	datFile = args.datFile
	oxFiles = args.oxFiles
	clusterFile = args.clusterFile
	nstep_skip = args.nstep_skip
	coarse_time = args.coarse_time
	center = args.center

	### check for conflicting inputs
	if geoFile is not None and oxFiles is not None:
		print("Error: Must provide either geometry file or oxDNA files, not both.")
		sys.exit()

	### output file parameters
	outFold = "lammpified/"
	outGeoFile = outFold + "geometry.in"
	outDatFile = outFold + "trajectory.dat"

	### create output folder
	ars.createEmptyFold(outFold)

	### lammpify lammps
	if geoFile is not None:

		### read geometry
		output = ars.readGeo(geoFile,getDbox3=True)
		if len(output) == 6:
			points_init,molecules,types,bonds,_,dbox3 = output
		elif len(output) == 7:
			points_init,molecules,types,_,bonds,_,dbox3 = output

		### read trajectory
		if datFile is not None:
			points,col2s,dbox3s = ars.readAtomDump(datFile,nstep_skip,coarse_time,getDbox3s=True)
			if center: points = ars.centerPointsMolecule(points,molecules,dbox3s,center='com',unwrap=False)

		### set colors
		if clusterFile is not None:
			npoint = points_init.shape[0]
			clusters = ars.readCluster(clusterFile)
			colors = getMoleculesFromClustersB1(clusters,npoint)
		else:
			colors = types

		### write output
		ars.writeGeo(outGeoFile,dbox3,points_init,types=colors,bonds=bonds,natomType=max(colors))
		if datFile is not None:
			ars.writeAtomDump(outDatFile,dbox3s,points,colors)

	### lammpify oxdna
	elif oxFiles is not None:

		### read topology (required)
		strands,bonds,nba_total = readTop(oxFiles[0])

		### read and center trajecotry
		points,dbox3 = ars.readOxDNA(oxFiles[1],nstep_skip,coarse_time,getDbox3=True)
		if center: points = ars.centerPointsMolecule(points,strands,dbox3,center='com',unwrap=True)

		### set colors
		if clusterFile is not None:
			clusters = ars.readCluster(clusterFile)
			colors = getMoleculesFromClustersB0(clusters,nba_total)
		else:
			strand_scaffold = stats.mode(strands).mode
			colors = np.where(strands == strand_scaffold, 1, 2)

		### write output
		dbox3s = np.ones((points.shape[0],3))*dbox3
		ars.writeGeo(outGeoFile,dbox3s[0],points[0,:,:],types=colors,bonds=bonds,natomType=max(colors))
		if points.shape[0] > 1:
			ars.writeAtomDump(outDatFile,dbox3s,points,colors)

	### error
	else:
		print("Error: Unknown simulation type, try again.")
		sys.exit()


################################################################################
### File Handlers

### read oxdna topology
def readTop(topFile):
	ars.checkFileExist(topFile,"topology")
	with open(topFile) as f:
		content = f.readlines()
	nba_total = int(content[0].split()[0])
	strands = np.zeros(nba_total,dtype=int)
	bonds = np.ones((nba_total,3),dtype=int)
	bond_count = 0
	for i in range(nba_total):
		strands[i] = int(content[i+1].split()[0])
		if int(content[i+1].split()[3]) != -1:
			bonds[bond_count,1] = i+1
			bonds[bond_count,2] = int(content[i+1].split()[3])+1
			bond_count += 1
	bonds = bonds[:bond_count,:]
	return strands, bonds, nba_total


################################################################################
### Utility Functions

### identify molecules (1 for unidentified, 2+ for molecule IDs)
def getMoleculesFromClustersB0(clusters,npoint):
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
def getMoleculesFromClustersB1(clusters,npoint):
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

