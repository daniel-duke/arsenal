import arsenal as ars
import numpy as np
import argparse
import sys

## Description
# this script reads an oxDNA topology (.top) and configuration (.dat) file,
  # identifies the scaffold (the longest strand), and rewrites both files with
  # the scaffold moved to the front, renumbered as strand 1.
# strand and base ordering are otherwise preserved: the scaffold is relocated
  # as a contiguous block at the start of the file, and the remaining strands
  # keep their relative order (renumbered 2, 3, ... in order of appearance).
# base-pair connectivity (3' and 5' neighbor indices) is remapped to match the
  # new base ordering.
# this script runs from the command line with two arguments: (1) the topology
  # file, (2) the configuration file.
# the script outputs files with the same names as the inputs, modified with
  # "_reorder" appended (but before the extension).


################################################################################
### Heart

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('topFile', type=str)
	parser.add_argument('datFile', type=str)
	args = parser.parse_args()

	### read topology and configuration
	nba_total, nstrand, strands, bases, n3s, n5s = readTop(args.topFile)
	points, dbox, axes = ars.readOxDNA(args.datFile, getAxes=True, keepOxUnits=True, ignorePBC=True)

	### check consistency
	if points.shape[1] != nba_total:
		print("Error: Topology and configuration files contain different numbers of bases.\n")
		sys.exit()

	### put scaffold first
	order, strands_new, bases_new, n3s_new, n5s_new = putScaffoldFirst(strands, bases, n3s, n5s)

	### write reordered topology
	outTopFile = ars.addSuffix(args.topFile, "_reorder")
	writeTop(outTopFile, nba_total, nstrand, strands_new, bases_new, n3s_new, n5s_new)

	### write reordered configuration
	outDatFile = ars.addSuffix(args.datFile, "_reorder")
	ars.writeOxDNA(outDatFile, points[:,order], axes[:,order], dbox)


################################################################################
### File Handlers

### read oxDNA topology file
def readTop(topFile):
	ars.checkFileExist(topFile, "topology")
	with open(topFile, 'r') as f:
		content = f.readlines()
	nba_total, nstrand = [ int(x) for x in content[0].split()[:2] ]
	strands = np.zeros(nba_total, dtype=int)
	bases = [None]*nba_total
	n3s = np.zeros(nba_total, dtype=int)
	n5s = np.zeros(nba_total, dtype=int)
	for i in range(nba_total):
		line = content[i+1].split()
		strands[i] = int(line[0])
		bases[i] = line[1]
		n3s[i] = int(line[2])
		n5s[i] = int(line[3])
	return nba_total, nstrand, strands, bases, n3s, n5s


### write oxDNA topology file
def writeTop(topFile, nba_total, nstrand, strands, bases, n3s, n5s):
	with open(topFile, 'w') as f:
		f.write(f"{nba_total} {nstrand}\n")
		for i in range(nba_total):
			f.write(f"{strands[i]} {bases[i]} {n3s[i]} {n5s[i]}\n")


################################################################################
### Utility Functions

### reorder bases so the longest strand (the scaffold) comes first, renumbered as strand 1
def putScaffoldFirst(strands, bases, n3s, n5s):
	nba_total = len(strands)

	### identify scaffold as the longest strand
	strand_ids, counts = np.unique(strands, return_counts=True)
	scaffold_id = strand_ids[np.argmax(counts)]
	scaffold_len = int(counts.max())

	### order bases: scaffold block first, then the rest in their original order
	is_scaffold = strands == scaffold_id
	order = np.concatenate([ np.flatnonzero(is_scaffold), np.flatnonzero(~is_scaffold) ])

	### map old base indices to new base indices
	old_to_new = np.zeros(nba_total, dtype=int)
	old_to_new[order] = np.arange(nba_total)

	### renumber strands: scaffold becomes 1, others keep relative order starting at 2
	other_ids = []
	for sid in strands[~is_scaffold]:
		if sid not in other_ids:
			other_ids.append(sid)
	strand_map = { scaffold_id: 1 }
	for i,sid in enumerate(other_ids):
		strand_map[sid] = i+2

	### build reordered arrays
	strands_new = np.array([ strand_map[s] for s in strands[order] ])
	bases_new = [ bases[i] for i in order ]
	n3s_new = np.array([ old_to_new[n] if n != -1 else -1 for n in n3s[order] ])
	n5s_new = np.array([ old_to_new[n] if n != -1 else -1 for n in n5s[order] ])

	### result
	return order, strands_new, bases_new, n3s_new, n5s_new


### run the script
if __name__ == "__main__":
	main()
	print()
