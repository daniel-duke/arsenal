import arsenal as ars
import numpy as np
from scipy import stats
import argparse
import sys

## Description
# this script reads an oxDNA topology (.top) file, two oxDNA configuration
  # (.dat) files, and a cluster file, and moves each cluster in the second
  # (non-ideal) configuration to match the scaffold positions in the first
  # (ideal) configuration.
# the scaffold must be the same in both configurations, and it must be the
  # first strand in both (the topology is used only to get the scaffold length);
  # the staples may differ between the two configurations.
# the cluster file belongs to the second configuration, with base indices
  # starting from 0; every nucleotide must belong to exactly one cluster, and
  # every cluster must contain at least 3 scaffold nucleotides.
# for each cluster, the scaffold nucleotides in the second configuration are
  # aligned to their ideal positions (Kabsch algorithm), and the resulting
  # rigid-body transformation is applied to the entire cluster (positions and
  # orientations).
# this script runs from the command line with four arguments: (1) the topology
  # file (for the second configuration), (2) the ideal configuration file, (3)
  # the non-ideal configuration file, (4) the cluster file.
# the script outputs a configuration file with the same name as the non-ideal
  # configuration file, modified with "_matched" appended (but before the
  # extension).


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('topFile', type=str)
	parser.add_argument('datFile_ideal', type=str)
	parser.add_argument('datFile_real', type=str)
	parser.add_argument('clusterFile', type=str)
	args = parser.parse_args()

	### read topology and get scaffold length
	strands, nba_total = readTop(args.topFile)
	nba_scaf = getScaffoldLength(strands)

	### read configurations
	points_ideal = ars.readOxDNA(args.datFile_ideal, keepOxUnits=True, ignorePBC=True)[0]
	points_real, dbox, axes_real = ars.readOxDNA(args.datFile_real, getAxes=True, keepOxUnits=True, ignorePBC=True)
	points_ideal = points_ideal[0]
	points_real = points_real[0]
	axes_real = axes_real[0]

	### check consistency
	if points_real.shape[0] != nba_total:
		print("Error: Topology and non-ideal configuration files contain different numbers of bases.\n")
		sys.exit()
	if points_ideal.shape[0] < nba_scaf:
		print("Error: Ideal configuration file contains fewer bases than the scaffold.\n")
		sys.exit()

	### read and check clusters
	clusters = ars.readCluster(args.clusterFile)
	checkClusters(clusters, nba_total, nba_scaf)

	### match clusters to ideal scaffold
	points_matched, axes_matched = matchClusters(points_real, axes_real, points_ideal, clusters, nba_scaf)

	### write matched configuration
	outDatFile = ars.addSuffix(args.datFile_real, "_matched")
	ars.writeOxDNA(outDatFile, points_matched, axes_matched, dbox)


################################################################################
### File Handlers

### read oxDNA topology file
def readTop(topFile):
	ars.checkFileExist(topFile, "topology")
	with open(topFile, 'r') as f:
		content = f.readlines()
	nba_total = int(content[0].split()[0])
	strands = np.zeros(nba_total, dtype=int)
	for i in range(nba_total):
		strands[i] = int(content[i+1].split()[0])
	return strands, nba_total


################################################################################
### Utility Functions

### get scaffold length, ensuring the scaffold is the first strand
def getScaffoldLength(strands):
	strand_scaffold = stats.mode(strands).mode
	if strand_scaffold != strands[0]:
		print("Error: Scaffold (most common strand) is not the first strand.\n")
		sys.exit()
	nba_scaf = int(np.sum(strands == strand_scaffold))
	if not np.all(strands[:nba_scaf] == strand_scaffold):
		print("Error: Scaffold nucleotides are not contiguous at the start of the topology.\n")
		sys.exit()
	return nba_scaf


### ensure every nucleotide is in exactly one cluster, and every cluster has enough scaffold
def checkClusters(clusters, nba_total, nba_scaf):

	### count cluster memberships
	counts = np.zeros(nba_total, dtype=int)
	for c in range(len(clusters)):
		for bai in clusters[c]:
			if bai < 0 or bai >= nba_total:
				print(f"Error: Cluster {c+1} contains base index {bai}, which is outside the origami (0 to {nba_total-1}).\n")
				sys.exit()
			counts[bai] += 1

	### check memberships
	bais_none = np.flatnonzero(counts == 0)
	bais_multi = np.flatnonzero(counts > 1)
	if len(bais_none) > 0:
		print(f"Error: {len(bais_none)} nucleotides are not in any cluster (first: {bais_none[0]}).\n")
		sys.exit()
	if len(bais_multi) > 0:
		print(f"Error: {len(bais_multi)} nucleotides are in multiple clusters (first: {bais_multi[0]}).\n")
		sys.exit()

	### check scaffold counts
	for c in range(len(clusters)):
		nscaf = sum(1 for bai in clusters[c] if bai < nba_scaf)
		if nscaf < 3:
			print(f"Error: Cluster {c+1} contains only {nscaf} scaffold nucleotides (at least 3 required for alignment).\n")
			sys.exit()


### align each cluster to the ideal scaffold positions
def matchClusters(points_real, axes_real, points_ideal, clusters, nba_scaf):
	points_matched = np.zeros_like(points_real)
	axes_matched = np.zeros_like(axes_real)
	for c in range(len(clusters)):

		### order cluster bases with scaffold first
		bais = np.array(clusters[c], dtype=int)
		bais_scaf = bais[bais < nba_scaf]
		bais = np.concatenate([ bais_scaf, bais[bais >= nba_scaf] ])

		### align cluster using scaffold bases
		r_aligned, R = ars.kabschAlgorithm(points_real[bais], points_ideal[bais_scaf], indices=np.arange(len(bais_scaf)), getR=True)
		points_matched[bais] = r_aligned
		axes_matched[bais] = axes_real[bais] @ R

	### result
	return points_matched, axes_matched


### run the script
if __name__ == "__main__":
	main()
	print()
