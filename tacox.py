import arsenal as ars
import numpy as np
import random
import argparse
import json
import sys


################################################################################
### Parameters

### physical constants
pitch = 0.332
d_com_axis = 0.6*ars.ox2nm

### start
def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('cadFile',type=str)
	parser.add_argument('lattice',type=str,nargs='?',default=None)
	parser.add_argument('--seed',type=int,default=42)
	parser.add_argument('--center',action='store_true')
	args = parser.parse_args()

	### determine lattice
	lattice = args.lattice
	if lattice is None:
		lattice = guessLattice(args.cadFile)

	### constants for square lattice
	if lattice == 'sq':
		class geo:
			lattice = 'sq'
			turns_per_domain = 3/32*8
			nnt_per_domain = 8
			d_vstrand = 2.6*ars.ox2nm

	### constants for hexagonal lattice
	elif lattice == 'hex':
		class geo:
			lattice = 'hex'
			turns_per_domain = 2/21*7
			nnt_per_domain = 7
			d_vstrand = 2.55*ars.ox2nm

	### error
	else:
		print("Error: Unrecognized lattice type.\n")
		sys.exit()


################################################################################
### Heart

	### set random seed
	random.seed(args.seed)

	### read caDNAno file
	topology, seq, occupancy_scaf, occupancy_stap, vstrand_ijs = parseCaDNAno(args.cadFile)

	### set positions
	r, axes, dbox = setPositions(occupancy_scaf, occupancy_stap, vstrand_ijs, geo)
	if args.center:	r = ars.centerPointsBead(r, dbox, report=False)

	### flip strands
	seq, r, axes = flipPrimeDir(topology, seq, r, axes)

	### write topology
	topFile = ars.changeExtension(args.cadFile, "top")
	writeTop(topFile, topology, seq)

	### write configuration
	print("Writing configuration...")
	datFile = ars.changeExtension(args.cadFile, "dat")
	ars.writeOxDNA(datFile, r, axes, dbox)


################################################################################
### File Handlers

### use length of vstrand to guess lattice type
def guessLattice(cadFile):
	ars.checkFileExist(cadFile, "caDNAno")
	nbox_vstrand = len(json.load(open(cadFile))["vstrands"][0]["scaf"])

	### square
	if nbox_vstrand % 32 == 0 and not nbox_vstrand % 21 == 0:
		print("Using square lattice.")
		return "sq"

	### hexagonal
	elif nbox_vstrand % 21 == 0 and not nbox_vstrand % 32 == 0:
		print("Using hexagonal lattice.")
		return "hex"

	### error
	else:
		print("Error: Cannot guess lattice type, please supply as second argument.")
		sys.exit()


### write topology file
def writeTop(topFile, topology, seq):
	print("Writing topology...")
	nba_total = len(topology)
	nstrand = max(nt[0] for nt in topology)+1
	with open(topFile, 'w') as f:
		f.write(f"{nba_total} {nstrand}\n")
		for i in range(nba_total):
			f.write(f"{topology[i][0]+1} {seq[i]} {topology[i][1]} {topology[i][2]}\n")


### extract necessary info from caDNAno file
def parseCaDNAno(cadFile):
	print("Parsing caDNAno file...")
	
	### load caDNAno file
	ars.checkFileExist(cadFile, "caDNAno")
	with open(cadFile, 'r') as f:
		json_string = f.read()
	j = json.loads(json_string)

	### get vstrand grid indices
	vstrand_ijs = parseVstrands(j)
	
	### get topology and the nucleotide(s) in each box
	topology_scaf, occupancy_scaf = parseStrands(j, "scaf")
	print(len(topology_scaf))
	topology_stap, occupancy_stap = parseStrands(j, "stap")

	### set scaffold sequence
	nnt_scaf = len(topology_scaf)
	seq_scaf = random.choices('ATCG', k=nnt_scaf)

	### set staple sequence
	nnt_stap = len(topology_stap)
	seq_stap = setStapSeq(seq_scaf, occupancy_scaf, occupancy_stap, nnt_stap)

	### combine scaffold and staple data
	nstrand_scaf = topology_scaf[-1][0]+1
	for nt in topology_stap:
		nt[0] += nstrand_scaf
		if nt[1] >=0: nt[1] += nnt_scaf
		if nt[2] >=0: nt[2] += nnt_scaf
	topology = topology_scaf + topology_stap
	seq = seq_scaf + seq_stap

	### result
	return topology, seq, occupancy_scaf, occupancy_stap, vstrand_ijs


### work function for parsing caDNAno
def parseStrands(j, strand_type):

	### count
	nvstrand = len(j["vstrands"])
	nbox_vstrand = len(j["vstrands"][0][strand_type])

	### initialize
	boxes_seen = set()
	si_curr = -1
	ni_curr = -1
	topology = []
	occupancy = [ [ [] for j in range(nbox_vstrand) ] for i in range(nvstrand) ]

	### loop over vstrands
	for vi in range(nvstrand):
		vstrand = j["vstrands"][vi]

		### loop over boxes
		for bi in range(nbox_vstrand):

			### check if nucleotide exists in this box
			if vstrand[strand_type][bi][0] == -1 and vstrand[strand_type][bi][2] == -1:
				continue

			### check if nucleotide has been logged
			box = vi*nbox_vstrand + bi
			if box in boxes_seen:
				continue

			### new strand found
			si_curr += 1

			### search for starting point
			vi_curr = vi
			bi_curr = bi
			while True:
				vi_next = j["vstrands"][vi_curr][strand_type][bi_curr][0]
				bi_next = j["vstrands"][vi_curr][strand_type][bi_curr][1]

				### found 5p end
				if vi_next == -1:
					vi_start = vi_curr
					bi_start = bi_curr
					break

				### back at original nucleotide (circular strand)
				elif vi_next == vi and bi_next == bi:
					vi_start = vi
					bi_start = bi
					break

				### keep searching
				else:
					vi_curr = vi_next
					bi_curr = bi_next

			### take note of starting nucleotide index (for circular strands)
			ni_start = ni_curr + 1

			### walk along strand (5p to 3p)
			start = True
			vi_curr = vi_start
			bi_curr = bi_start
			while True:

				### log box
				boxes_seen.add(vi_curr*nbox_vstrand + bi_curr)

				### log skip
				if j["vstrands"][vi_curr]["skip"][bi_curr] != 0:
					occupancy[vi_curr][bi_curr].append(-1)

				### new nucleotide
				else:
					ni_curr += 1
					occupancy[vi_curr][bi_curr].append(ni_curr)
					topology.append([si_curr,ni_curr-1,ni_curr+1])

					### assume 5p end empty for now (fixed later for circular strands)
					if start:
						topology[ni_curr][1] = -1
						start = False

					### loop adjustments
					for li in range(j["vstrands"][vi_curr]["loop"][bi_curr]):
						ni_curr += 1
						occupancy[vi_curr][bi_curr].append(ni_curr)
						topology.append([si_curr,ni_curr-1,ni_curr+1])

				### peek at next nucleotide
				vi_next = j["vstrands"][vi_curr][strand_type][bi_curr][2]
				bi_next = j["vstrands"][vi_curr][strand_type][bi_curr][3]

				### found 3p end
				if vi_next == -1:
					topology[ni_curr][2] = -1
					break

				### back at start (circular strand)
				elif vi_next == vi_start and bi_next == bi_start:
					topology[ni_curr][2] = ni_start
					topology[ni_start][1] = ni_curr
					break

				### keep on going
				else:
					vi_curr = vi_next
					bi_curr = bi_next

	### result
	return topology, occupancy


### get vstrand grid location
def parseVstrands(j):
	nvstrand = len(j["vstrands"])
	vstrand_ijs = np.zeros((nvstrand,2))
	for vi in range(nvstrand):
		vstrand_ijs[vi,0] = j["vstrands"][vi]["col"]
		vstrand_ijs[vi,1] = j["vstrands"][vi]["row"]
	return vstrand_ijs


################################################################################
### Positions

### calculate positions of all nucleotides
def setPositions(occupancy_scaf, occupancy_stap, vstrand_ijs, geo):
	print("Setting positions...")

	### count
	nnt_scaf = max(x for vstrand in occupancy_scaf for box in vstrand for x in box)+1
	nnt_stap = max(x for vstrand in occupancy_stap for box in vstrand for x in box)+1

	### initialize
	r_scaf = np.zeros((nnt_scaf,3))
	r_stap = np.zeros((nnt_stap,3))
	axes_scaf = np.zeros((nnt_scaf,3,3))
	axes_stap = np.zeros((nnt_stap,3,3))

	### count
	nvstrand = len(occupancy_scaf)
	ndom_vstrand = len(occupancy_scaf[0])//geo.nnt_per_domain

	### loop over virtual strands
	for vi in range(nvstrand):

		### get scaf direction
		scaf5to3 = isScaf5to3(vstrand_ijs[vi])

		### loop over domains
		for di in range(ndom_vstrand):

			### figure out nucleotide placement within domain
			placement, nplace = mapBoxToPlacement(occupancy_scaf, occupancy_stap, vi, di, geo)

			### no nucleotides to place
			if all(not x for x in placement):
				continue

			### get positions for each place
			r_domain_scaf, r_domain_stap, axes_domain_scaf, axes_domain_stap = getPositionsDomain(nplace, di, vstrand_ijs[vi], geo)

			### loop over boxes
			for bi_domain in range(geo.nnt_per_domain):
				bi = di*geo.nnt_per_domain + bi_domain
				nnt_box = len(placement[bi_domain])

				### loop over nucleotides in box
				for ni_box in range(nnt_box):

					### determine direction of placement within box
					if scaf5to3:
						ni_box_scaf = ni_box
						ni_box_stap = nnt_box-ni_box-1
					else:
						ni_box_scaf = nnt_box-ni_box-1
						ni_box_stap = ni_box

					### set scaffold positions
					if occupancy_scaf[vi][bi]:
						ni = occupancy_scaf[vi][bi][ni_box_scaf]
						r_scaf[ni] = r_domain_scaf[placement[bi_domain][ni_box]]
						axes_scaf[ni] = axes_domain_scaf[placement[bi_domain][ni_box]]

					### set staple positions
					if occupancy_stap[vi][bi]:
						ni = occupancy_stap[vi][bi][ni_box_stap]
						r_stap[ni] = r_domain_stap[placement[bi_domain][ni_box]]
						axes_stap[ni] = axes_domain_stap[placement[bi_domain][ni_box]]

	### convert to oxDNA units
	r_scaf *= ars.nm2ox
	r_stap *= ars.nm2ox

	### combine scaffold and staple
	r = np.concatenate((r_scaf, r_stap), axis=0)
	axes = np.concatenate((axes_scaf, axes_stap), axis=0)

	### determine box size
	dbox = 2*max( np.max(r,axis=0) - np.min(r,axis=0) )

	### result
	return r, axes, dbox


### figure out spacing for nucleotides within domain
def mapBoxToPlacement(occupancy_scaf, occupancy_stap, vi, di, geo):

	### full domain (excluding skips, including loops)
	domain = []

	### loop over boxes
	bi_start = di*geo.nnt_per_domain
	for bi in range(bi_start,bi_start+geo.nnt_per_domain):

		### count
		nnt_box = max([ len(occupancy_scaf[vi][bi]), len(occupancy_stap[vi][bi]) ])

		### empty
		if nnt_box == 0:
			domain.append(None)

		### skip
		elif occupancy_scaf[vi][bi] and occupancy_scaf[vi][bi][0] == -1:
			continue

		### nucleotides
		else:
			for ni_box in range(nnt_box):
				domain.append([bi-bi_start,ni_box])

	### remove empty places if necessary and possible
	i = 0
	while len(domain) > geo.nnt_per_domain and i < len(domain):
		if domain[i] is None:
			domain.pop(i)
		else:
			i += 1

	### set placement
	nplace = len(domain)
	placement = [ [] for i in range(geo.nnt_per_domain) ]
	for pi in range(len(domain)):
		if domain[pi] == None:
			continue
		bi_domain = domain[pi][0]
		placement[bi_domain].append(pi)

	### result
	return placement, nplace


### calculate positions for nucleotides in domain
def getPositionsDomain(nplace, di, vstrand_ij, geo):

	### initialize
	r_scaf = np.zeros((nplace,3))
	r_stap = np.zeros((nplace,3))
	axes_scaf = np.zeros((nplace,3,3))
	axes_stap = np.zeros((nplace,3,3))

	### determine prime direction
	scaf5to3 = isScaf5to3(vstrand_ij)

	### calculate starting angle (xy)
	if geo.lattice == 'sq':
		if scaf5to3:
			theta_start = 178 + di*geo.turns_per_domain*360
		else:
			theta_start = -2 + di*geo.turns_per_domain*360
	else:
		if scaf5to3:
			theta_start = -37 + di*geo.turns_per_domain*360
		else:
			theta_start = 143 + di*geo.turns_per_domain*360

	### calculate starting z-axis position
	z_start = di*geo.nnt_per_domain*pitch

	### loop over places
	for i in range(nplace):

		### get angle
		theta = theta_start + (i+0.5)*geo.turns_per_domain/nplace*360
		theta_rad = np.deg2rad(theta)

		### get rotation matrix
		R = np.zeros((2,2))
		R[0,0] = np.cos(theta_rad)
		R[1,0] = -np.sin(theta_rad)
		R[0,1] = np.sin(theta_rad)
		R[1,1] = np.cos(theta_rad)

		### set scaffold axes
		axes_scaf[i] = np.eye(3)
		if scaf5to3:
			axes_scaf[i,1:] *= -1
		axes_scaf[i,0,:2] = -axes_scaf[i,0,:2] @ R
		axes_scaf[i,1,:2] = -axes_scaf[i,1,:2] @ R

		### set staple axes
		axes_stap[i] = np.eye(3)
		if not scaf5to3:
			axes_stap[i,1:] *= -1
		axes_stap[i,0,:2] = axes_stap[i,0,:2] @ R
		axes_stap[i,1,:2] = axes_stap[i,1,:2] @ R

		### get helix xy position
		if geo.lattice == "sq":
			vstrand_xy = geo.d_vstrand*vstrand_ij
		else:
			vstrand_xy = np.zeros(2)
			vstrand_xy[0] = vstrand_ij[0]*geo.d_vstrand*np.sqrt(3)/2
			vstrand_xy[1] = vstrand_ij[1]*geo.d_vstrand*1.5
			if not scaf5to3:
				vstrand_xy[1] += geo.d_vstrand*0.5

		### set xy position
		r_scaf[i,:2] = vstrand_xy - d_com_axis*axes_scaf[i,0,:2]
		r_stap[i,:2] = vstrand_xy - d_com_axis*axes_stap[i,0,:2]

		### set z-axis positions
		r_scaf[i,2] = z_start + i*pitch*geo.nnt_per_domain/nplace
		r_stap[i,2] = z_start + i*pitch*geo.nnt_per_domain/nplace

	### result
	return r_scaf, r_stap, axes_scaf, axes_stap


### determine prime direction
def isScaf5to3(vstrand_ij):
	return int((sum(vstrand_ij)+1)%2)


################################################################################
### Sequence

### flip strands (5to3 -> 3to5, and vice versa)
def flipPrimeDir(topology, seq, r, axes):

	### count
	nnt_total = len(seq)
	nstrand = max(nt[0] for nt in topology)+1

	### initialize
	seq_flip = [None]*nnt_total
	r_flip = np.zeros((nnt_total,3))
	axes_flip = np.zeros((nnt_total,3,3))

	### loop over strands
	ni_start = 0
	for si in range(nstrand):
		nnt_strand = sum(1 for nt in topology if nt[0]==si)

		### loop over nucleotides
		for ni_strand in range(nnt_strand):

			### calculate new nucleotide
			ni = ni_start+ni_strand
			ni_flip = ni_start+nnt_strand-ni_strand-1

			### flippy flip flip
			seq_flip[ni_flip] = seq[ni]
			r_flip[ni_flip] = r[ni]
			axes_flip[ni_flip] = axes[ni]

		### next strand starting point
		ni_start += nnt_strand

	### result
	return seq_flip, r_flip, axes_flip


### make staples complementary to scaffold, random if no complement
def setStapSeq(seq_scaf, occupancy_scaf, occupancy_stap, nnt_stap):
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

	### count
	nvstrand = len(occupancy_stap)
	nbox = len(occupancy_stap[0])

	### loop over vstrands
	seq_stap = [None]*nnt_stap
	for vi in range(nvstrand):

		### loop over boxes
		for bi in range(nbox):
			nnt_box = len(occupancy_stap[vi][bi])

			### loop over nucleotides within box
			for ni_box in range(nnt_box):
				ni_stap = occupancy_stap[vi][bi][ni_box]

				### set as complement
				if occupancy_scaf[vi][bi]:
					ni_box_rev = nnt_box-ni_box-1
					ni_scaf = occupancy_scaf[vi][bi][ni_box_rev]
					seq_stap[ni_stap] = comp[seq_scaf[ni_scaf]]

				### set randomly
				else:
					seq_stap[ni_stap] = random.choice('ATCG')

	### result				
	return seq_stap


### run the script
if __name__ == "__main__":
	main()
	print()


