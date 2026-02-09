import arsenal as ars
import numpy as np
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
	parser.add_argument('--sim-type',		type=str, required=True,	help="simulation type (lammps, oxdna)")
	parser.add_argument('--geoFile',		type=str, default=None,		help="name of lammps geometry file")
	parser.add_argument('--topFile',		type=str, default=None,		help="name of oxdna topology file")
	parser.add_argument('--datFile',		type=str, default=None,		help="name of trajectory file")
	parser.add_argument('--clusterFile',	type=str, default=None,		help="name of clusters file")
	parser.add_argument('--nstep-skip',		type=int, default=0, 		help="number of recorded initial steps to skip")
	parser.add_argument('--coarse-time',	type=int, default=1, 		help="coarse factor for time steps")
	parser.add_argument('--center',			action='store_true',		help="whether to center the trajectory")

	### set arguments
	args = parser.parse_args()
	sim_type = args.sim_type
	geoFile = args.geoFile
	topFile = args.topFile
	datFile = args.datFile
	clusterFile = args.clusterFile
	nstep_skip = args.nstep_skip
	coarse_time = args.coarse_time
	center = args.center

	### output file parameters
	outFold = "lammpified/"
	outGeoFile = outFold + "geometry.in"
	outDatFile = outFold + "trajectory.dat"

	### create output folder
	ars.createEmptyFold(outFold)

	### lammpify lammps
	if sim_type == "lammps":

		### read geometry (required)
		if geoFile is not None:
			output = ars.readGeo(geoFile,getDbox3=True)
			if len(output) == 6:
				points_init,molecules,types,bonds,_,dbox3 = output
			elif len(output) == 7:
				points_init,molecules,types,_,bonds,_,dbox3 = output
		else:
			print("Error: geometry file required for lammps simulations.\n"); sys.exit()

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
			writeAtomDump(outDatFile,dbox3s,points,colors)

	### lammpify oxdna
	elif sim_type == "oxdna":

		### read topology (required)
		if topFile is not None:
			strands,bonds,nba_total = readTop(topFile)
		else:
			print("Error: topology file required for oxdna simulations.\n"); sys.exit()

		### read and center trajecotry (required)
		if datFile is not None:
			points,dbox3 = ars.readOxDNA(datFile,nstep_skip,coarse_time,getDbox3=True)
			if center: points = ars.centerPointsMolecule(points,strands,dbox3,center='com',unwrap=True)
		else:
			print("Error: trajectory file required for oxdna simulations.\n"); sys.exit()

		### set colors
		if clusterFile is not None:
			clusters = ars.readCluster(clusterFile)
			colors = getMoleculesFromClustersB0(clusters,nba_total)
		else:
			colors = np.minimum(strands,2)

		### write output
		dbox3s = np.ones((points.shape[0],3))*dbox3
		ars.writeGeo(outGeoFile,dbox3s[0],points[0,:,:],types=colors,bonds=bonds,natomType=max(colors))
		if points.shape[0] > 1:
			writeAtomDump(outDatFile,dbox3s,points,colors)

	### error
	else:
		print("Error: unknown simulation type, try again.")
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
	return strands,bonds,nba_total


### write lammps-style atom dump
def writeAtomDump(outDatFile,dbox3s,points,colors):
	print("Writing trajectory...")
	nstep = points.shape[0]
	npoint = points.shape[1]
	len_npoint = len(str(npoint))
	len_ncolor = len(str(max(colors)))
	with open(outDatFile,'w') as f:
		for i in range(nstep):
			len_dbox = len(str(int(max(dbox3s[i]))))
			f.write(f"ITEM: TIMESTEP\n{i}\n")
			f.write(f"ITEM: NUMBER OF ATOMS\n{npoint}\n")
			f.write(f"ITEM: BOX BOUNDS pp pp pp\n")
			f.write(f"-{dbox3s[i,0]/2:0{len_dbox+3}.2f} {dbox3s[i,0]/2:0{len_dbox+3}.2f} xlo xhi\n")
			f.write(f"-{dbox3s[i,1]/2:0{len_dbox+3}.2f} {dbox3s[i,1]/2:0{len_dbox+3}.2f} ylo yhi\n")
			f.write(f"-{dbox3s[i,2]/2:0{len_dbox+3}.2f} {dbox3s[i,2]/2:0{len_dbox+3}.2f} zlo zhi\n")
			f.write("ITEM: ATOMS id type xs ys zs\n")
			for j in range(npoint):
				f.write(f"{j+1:<{len_npoint}} " + \
						f"{colors[j]:<{len_ncolor}}  " + \
						f"{points[i,j,0]/dbox3s[i,0]+1/2:10.8f} " + \
						f"{points[i,j,1]/dbox3s[i,1]+1/2:10.8f} " + \
						f"{points[i,j,2]/dbox3s[i,2]+1/2:10.8f}\n")


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

