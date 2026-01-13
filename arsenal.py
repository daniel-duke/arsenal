import arsenal as ars
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import gaussian_kde
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.special import gammaln
import subprocess
import shutil
import pickle
import os
import sys
import re
import ast

## Description
# Welcome! You have stumbled upon my humble arsenal of python functions for
  # analyzing simulations relevant to DNA origami, from oxDNA to DNAfold to
  # mesoscopic LAMMPS models. Enjoy!

## Notes
# For all plotting functions, 

## Terminology
# bais - 2D python array containing lists of base indices (starting from 0)
# bdis - 2D python array containing lists of bead indices (starting from 1)
# r - npoint x 3 numpy array containing positions
# points - nstep x npoint x 3 numpy array containing positions
# molecules - npoint-element numpy array containing the molecule index for each point (starting from 1)
# col2s - npoint-element numpy array containing atomDump second column data for each point
# dbox - float, periodic box diameter
# dbox3 - 3-element array, periodic box diameter in each coordinate direction
# dbox3s - nstep x 3 array, periodic box diameter in each coordinate direction at each time step


################################################################################
### File Readers

### read oxdna trajectory
def readOxDNA(datFile, nstep_skip=0, coarse_time=1, bais='all', coarse_points=1, nstep_max='all', **kwargs):

	### additional keyword args
	checkBackbone 	= False 	if 'checkBackbone' not in kwargs else kwargs['checkBackbone']
	ntFirstStrand 	= 'auto'	if 'ntFirstStrand' not in kwargs else kwargs['ntFirstStrand']
	useHbondSite 	= False		if 'useHbondSite' not in kwargs else kwargs['useHbondSite']
	ignorePBC 		= False		if 'ignorePBC' not in kwargs else kwargs['ignorePBC']
	getDbox3 		= False		if 'useDbox3' not in kwargs else kwargs['useDbox3']
	getUsedEvery	= False		if 'getUsedEvery' not in kwargs else kwargs['getUsedEvery']
	report_every	= 1000		if 'report_every' not in kwargs else kwargs['report_every']

	### notes
	# assumes the bais array stores the base indices starting from 0.
	# assumes the the dump frequency does not change, which it never should.
	# assumes the box diameter does not change, which it never should.
	# unless instructed otherwise, returns only first dimension of box diameter.
	# all three dimensions of the box are used to get the points.
	# backbone bond checking assumes the the first ntFirstStrand nucleotides form a continous strand.
	# backbond bond checking assumes all bases form one continuous strand if ntFirstStrand is set to auto.

	### load trajectory file
	print("Loading oxDNA trajectory...")
	ars.checkFileExist(datFile, "trajectory")
	with open(datFile, 'r') as f:
		content = f.readlines()
	print("Parsing trajectory...")

	### extract metadata
	nba_total = 0
	while nba_total+3 < len(content) and ars.isnumber(content[nba_total+3].split()[0]):
		nba_total += 1
	dbox3 = np.array(content[1].split()[2:5],dtype=float)
	dump_every = 0
	if len(content) > nba_total+3:
		dump_every = int(content[nba_total+3].split()[2]) - int(content[0].split()[2])
	nstep_recorded = int(len(content)/(nba_total+3))
	nstep_trimmed = int((nstep_recorded-nstep_skip)/coarse_time)
	if nstep_trimmed <= 0:
		print("Error: Cannot read oxDNA trajectory - too much initial time cut off.\n")
		sys.exit()

	### determine number of steps to use
	if isinstance(nstep_max, str) and nstep_max == 'all':
		nstep_used = nstep_trimmed
	elif ars.isinteger(nstep_max):
		nstep_used = min([nstep_max,nstep_trimmed])
	else:
		print("Error: Cannot read oxDNA trajectory - max number of steps must be 'all' or integer.\n")
		sys.exit()

	### report step counts
	print("{:1.2e} steps in simulation".format(nstep_recorded*dump_every))
	print("{:1.2e} steps recorded".format(nstep_recorded))
	print("{:1.2e} steps used".format(nstep_used))

	### determine bases to use
	if isinstance(bais, str) and bais == 'all':
		bais = [range(int(np.ceil(nba_total/coarse_points)))]
	elif ars.isinteger(bais) and bais < 0:
		bdis = [range(int(np.ceil(-bais/coarse_points)))]
	elif ars.isinteger(bais):
		bais = [[bais]]
	elif ars.isarray(bais) and ars.isinteger(bais[0]):
		bais = [bais[::coarse_points]]
	elif ars.isarray(bais) and ars.isarray(bais[0]) and ars.isinteger(bais[0][0]):
		for i in range(len(bais)):
			bais[i] = bais[i][::coarse_points]
	else:
		print("Error: Cannot read oxDNA trajectory - base indices must be 'all', integer, 1D int array, or 2D int array.\n")
		sys.exit()

	### prepare for backbone bond checking
	if checkBackbone:
		if isinstance(ntFirstStrand, str) and ntFirstStrand == 'auto':
			ntFirstStrand = nba_total
		elif not ars.isinteger(ntFirstSrand):
			print("Error: Cannot read oxDNA trajectory - number of nucleotides in first strand must be 'auto' or integer.\n")
			sys.exit()

	### count total number of bases to use
	nba_used = 0
	for i in range(len(bais)):
		nba_used += len(bais[i])

	### extract the data
	points = np.zeros((nstep_used,nba_used,3))
	groups = np.zeros(nba_used, dtype=int)
	for i in range(nstep_used):
		point_count = 0
		for g in range(len(bais)):
			for j in range(len(bais[g])):
				bai = bais[g][j]
				if bai >= nba_total:
					print(f"Error: Cannot read oxDNA trajectory - requested base index {bai} exceeds the number of bases in the simulation ({nba_total}).\n")
					sys.exit()
				line = content[(nba_total+3)*(nstep_skip+i*coarse_time)+3+bai].split()
				if i == 0:
					groups[point_count] = g + 1
				if not useHbondSite:
					points[i,point_count] = np.array(line[:3],dtype=float)
				else:
					com = np.array(line[:3],dtype=float)
					a1 = np.array(line[3:6],dtype=float)
					points[i,point_count] = com + 0.4*a1
				point_count += 1
		if not ignorePBC:
			points[i] = ars.applyPBC( points[i], dbox3 )
		if (i+1)%report_every == 0:
			print(f"processed {i+1} steps...")

		### check backbone bond lengths
		if checkBackbone:
			r12_eq_max = 0.7564 + 0.25		# for oxDNA2, r_0_backbone + Delta_FENE
			stretched_count = 0
			for bai in range(ntFirstStrand-1):
				line_index = (nba_total+3)*(nstep_skip+i*coarse_time)+3+bai
				line0 = content[line_index+0].split()
				line1 = content[line_index+1].split()
				com_bai0 = np.array(line0[:3],dtype=float)
				com_bai1 = np.array(line1[:3],dtype=float)
				a1_bai0  = np.array(line0[3:6],dtype=float)
				a1_bai1  = np.array(line1[3:6],dtype=float)
				a3_bai0  = np.array(line0[6:9],dtype=float)
				a3_bai1  = np.array(line1[6:9],dtype=float)
				a2_bai0 = np.cross(a3_bai0, a1_bai0)
				a2_bai1 = np.cross(a3_bai1, a1_bai1)
				bb_bai0 = com_bai0 - 0.34*a1_bai0 + 0.3408*a2_bai0		# for oxDNA2, from old documentation
				bb_bai1 = com_bai1 - 0.34*a1_bai1 + 0.3408*a2_bai1		# for oxDNA2, from old documentation
				r12_eq  = np.linalg.norm(bb_bai1-bb_bai0)
				if r12_eq > r12_eq_max:
					stretched_count += 1

	### stretched backbone warning
	if checkBackbone and stretched_count > 0:
		print(f"Flag: Backbone bond exceeded max value {stretched_count} times.")
	
	### oxDNA units to nm
	points *= 0.8518
	dbox3 *= 0.8518

	### result
	output = [ points, dbox3[0] ]
	if getDbox3: output[-1] = dbox3
	if len(bais) > 1: output.append(groups)
	if getUsedEvery: output.append(dump_every*coarse_time)
	return output


### read lammps-style trajectory
def readAtomDump(datFile, nstep_skip=0, coarse_time=1, bdis='all', coarse_points=1, nstep_max='all', **kwargs):
	
	### additional keyword args
	ignorePBC		= False		if 'ignorePBC' not in kwargs else kwargs['ignorePBC']
	getDbox3s		= False		if 'getDbox3s' not in kwargs else kwargs['getDbox3s']
	getUsedEvery	= False		if 'getUsedEvery' not in kwargs else kwargs['getUsedEvery']
	report_every	= 1000		if 'report_every' not in kwargs else kwargs['report_every']
	
	### notes
	# assumes the bdis array stores the atom indices starting from 1.
	# assumes the dump frequency never changes, which it never should.
	# assumes there are five columns (id, col2, xs, ys, zs).
	# unless instructed otherwise, returns only first dimension of box diameter in first timestep.
	# all three dimensions of the box are extracted each timestep and used to get the points.
	# the returned points are centered about the origin, even if of the extracted values are not.

	### load trajectory file
	print("Loading LAMMPS-style trajectory...")
	ars.checkFileExist(datFile, "trajectory")
	with open(datFile, 'r') as f:
		content = f.readlines()
	print("Parsing trajectory...")

	### extract metadata
	nbd_total = int(content[3].split()[0])
	dump_every = 0
	if len(content) > nbd_total+10:
		dump_every = int(content[nbd_total+10].split()[0]) - int(content[1].split()[0])
	nstep_recorded = int(len(content)/(nbd_total+9))
	nstep_trimmed = int((nstep_recorded-nstep_skip-1)/coarse_time)+1
	if nstep_trimmed <= 0:
		print("Error: Cannot read atom dump - too much initial time cut off.\n")
		sys.exit()

	### interpret input
	if isinstance(nstep_max, str) and nstep_max == 'all':
		nstep_used = nstep_trimmed
	elif ars.isinteger(nstep_max):
		nstep_used = min([nstep_max,nstep_trimmed])
	else:
		print("Error: Cannot read atom dump - max number of steps must be 'all' or integer.\n")
		sys.exit()

	### report step counts
	print("{:1.2e} steps in simulation".format(nstep_recorded*dump_every))
	print("{:1.2e} steps recorded".format(nstep_recorded))
	print("{:1.2e} steps used".format(nstep_used))

	### interpret input
	if isinstance(bdis, str) and bdis == 'all':
		bdis = [list(range(1,int(np.ceil(nbd_total/coarse_points))+1))]
	elif ars.isinteger(bdis) and bdis < 0:
		bdis = [list(range(1,int(np.ceil(-bdis/coarse_points))+1))]
	elif ars.isinteger(bdis):
		bdis = [[bdis]]
	elif ars.isarray(bdis) and ars.isinteger(bdis[0]):
		bdis = [bdis[::coare_points]]
	elif ars.isarray(bdis) and ars.isarray(bdis[0]) and ars.isinteger(bdis[0][0]):
		for i in range(len(bdis)):
			bdis[i] = bdis[i][::coarse_points]
	else:
		print("Error: Cannot read atom dump - bead indices must be 'all', integer, 1D int array, or 2D int array.\n")
		sys.exit()

	### count total number of beads to use
	nbd_used = 0
	for i in range(len(bdis)):
		nbd_used += len(bdis[i])

	### extract the data
	points = np.zeros((nstep_used,nbd_used,3))
	col2s = np.zeros(nbd_used, dtype=int)
	groups = np.zeros(nbd_used, dtype=int)
	dbox3s = np.zeros((nstep_used,3))
	for i in range(nstep_used):
		point_count = 0
		for k in range(3):
			line = content[(nbd_total+9)*(nstep_skip+i*coarse_time)+5+k].split()
			dbox3s[i,k] = float(line[1]) - float(line[0])
		for g in range(len(bdis)):
			for j in range(len(bdis[g])):
				bdi = bdis[g][j]
				if bdi > nbd_total:
					print(f"Error: Cannot read atom dump - requested bead index {bdi} exceeds the number of beads in the simulation ({nbd_total}).\n")
					sys.exit()
				line = content[(nbd_total+9)*(nstep_skip+i*coarse_time)+9+bdi-1].split()
				if i == 0:
					col2s[point_count] = int(line[1])
					groups[point_count] = g + 1
				points[i,point_count] = (np.array(line[2:5],dtype=float)-1/2)*dbox3s[i]
				point_count += 1
		if not ignorePBC:
			points[i] = ars.applyPBC(points[i], dbox3s[i])
		if (i+1)%report_every == 0:
			print(f"processed {i+1} steps...")

	### results
	output = [ points, col2s ]
	output.append(dbox3s[0,0])
	if getDbox3s: output[-1] = dbox3s
	if len(bdis) > 1: output.append(groups)
	if getUsedEvery: output.append(dump_every*coarse_time)
	return output


### extract number of steps from lammps-style trajectory
def getNstep(datFile, nstep_skip=0, coarse_time=1):

	### load trajectory file
	ars.checkFileExist(datFile, "trajectory")
	with open(datFile, 'r') as f:
		content = f.readlines()

	### extract dump frequency
	nbd_total = int(content[3].split()[0])
	nstep_recorded = int(len(content)/(nbd_total+9))
	nstep_trimmed = int((nstep_recorded-nstep_skip-1)/coarse_time)+1
	return nstep_trimmed


### read the given columns from a dump file
def readDump(dumpFile, nstep_skip=0, coarse_time=1, cols='all', nstep_max='all', **kwargs):

	### additional keyword arguments
	isInt = False if 'isInt' not in kwargs else kwargs['isInt']

	### find dump file
	print("Parsing LAMMPS-style dump...")
	ars.checkFileExist(dumpFile, "dump")

	### get metadata
	with open(dumpFile, 'r') as f:
		content = iter(f)
		for line in range(9):
			next(content)
		ncol_dump = len(next(content).split())

	### interpret input
	if cols == 'all':
		cols = np.arange(ncol_dump,dtype=int)
	elif ars.isinteger(cols) and cols < ncol_dump:
		cols = np.array([cols])
	elif ars.isarray(cols) and ars.isinteger(cols[0]) and max(cols) < ncol_dump:
		cols = np.array(cols)
	else:
		print("Error: Cannot read dump - columns must be 'all', integer, or 1D integer array.")
		sys.exit()
	if isinstance(nstep_max, str) and nstep_max == 'all':
		nstep_max = 0
	elif not ars.isinteger(nstep_max):
		print("Error: Cannot read dump - max number of steps must be 'all' or integer.\n")
		sys.exit()

	### count columns
	ncol = len(cols)

	### initialize
	data = []
	step_count = 0

	### loop over lines
	with open(dumpFile, 'r') as f:
		content = iter(f)
		for line in content:

			### do nothing with first line and skip next 2 lines
			for line in range(2): next(content)

			### count entries for this step
			nentry = int(next(content))

			### skip next 5 lines
			for line in range(5): next(content)

			### determine if skipping step
			if step_count < nstep_skip or (step_count-nstep_skip)%coarse_time != 0:
				for line in range(nentry):
					next(content)

			### if reading step
			else:

				### initialize step data
				dtype = int if isInt else float
				data_step = np.zeros((nentry,ncol), dtype=dtype)

				### read data
				for i in range(nentry):
					line = next(content).split()
					data_step[i] = [line[c] for c in cols]

				### save step data
				data.append(data_step)

				### stop collecting data
				if len(data) == nstep_max:
					break

			### increase step count
			step_count = step_count + 1

	### result
	return data


### extract and plot thermo data from lammps output
def plotThermo(reportFile="report.out", cols=["E_mol","TotEng"]):

	### load report file
	ars.checkFileExist(reportFile, "report")
	with open(reportFile, 'r') as f:
		content = f.readlines()

	### find thermo headers
	idx_headers = [i for i, line in enumerate(content) if re.match(r'^\s*Step', line)]
	if not idx_headers:
		print("Error: Could not read thermo file - no thermo header (\"Step...\") found.")
	idx_header = idx_headers[-1]
	header = content[idx_header].split()

	### read the data
	data = []
	for line in content[idx_header+1:]:
		vals = line.split()
		if len(vals) != len(header) or not ars.isnumber(vals[0]):
			break
		vals = list(map(float,vals))
		data.append(vals)

	### process data
	if not data:
		print("Error: Could not read thermo file - no data found after header.")
		sys.exit()
	data = np.array(data)
	data_dict = {col: data[:,i] for i,col in enumerate(header)}
	for col in cols:
		if col not in data_dict:
			print(f"Warning: '{col}' not found in thermo columns, skipping.")
			continue

	### plot all data
	for col in cols:
		ars.plotLine(data_dict["Step"], data_dict[col], figLabel="Thermo", label=col)
	plt.xlabel("Time Step")
	plt.legend()

	### plot distributions
	for col in cols:
		ars.plotHist(data_dict[col], Alabel=col, figLabel=col, plotAvgLine=True, plotMedLine=True)


### read lammps-style geometry
def readGeo(geoFile, **kwargs):

	### keyword arguments
	extraLabel	= None	if 'extraLabel' not in kwargs else kwargs['extraLabel']
	getDbox		= None	if 'getDbox' not in kwargs else kwargs['getDbox']
	getDbox3	= None	if 'getDbox3' not in kwargs else kwargs['getDbox3']

	### notes
	# all indices are stored directly as they appear in the geometry file; in
	  # other words, the molecule indices, type indices, and atom indices in 
	  # the bonds all start from 1.

	### load geometry file
	ars.checkFileExist(geoFile, "geometry")
	with open(geoFile, 'r') as f:
		content = f.readlines()

	### get atom and bond counts
	natom = -1
	nbond = -1
	nangle = -1
	for i in range(len(content)):
		if len(content[i].split()) == 2:
			if content[i].split()[1] == 'atoms':
				natom = int(content[i].split()[0])
			if content[i].split()[1] == 'bonds':
				nbond = int(content[i].split()[0])
			if content[i].split()[1] == 'angles':
				nangle = int(content[i].split()[0])
		if natom != -1 and nbond != -1 and nangle != -1:
			break
		if i == len(content)-1:
			if natom == -1:
				natom = 0
			if nbond == -1:
				nbond = 0
			if nangle == -1:
				nangle = 0

	### get box diameter
	if getDbox or getDbox3:
		dbox3 = np.zeros(3)
		for i in range(len(content)):
			line = content[i].split()
			if len(line) == 4 and line[2] == 'xlo':
				dbox3[0] = float(line[1]) - float(line[0])
			if len(line) == 4 and line[2] == 'ylo':
				dbox3[1] = float(line[1]) - float(line[0])
			if len(line) == 4 and line[2] == 'zlo':
				dbox3[2] = float(line[1]) - float(line[0])
			if all(dbox3>0):
				break

	### get atom information
	r = np.zeros((natom,3))
	ids = np.zeros(natom, dtype=int)
	molecules = np.zeros(natom, dtype=int)
	types = np.zeros(natom, dtype=int)
	charges = np.zeros(natom)
	readCharge = False
	if natom:
		for i in range(len(content)):
			if len(content[i].split()) > 0 and content[i].split()[0] == 'Atoms':
				line_index = i+2
				break	
		for i in range(natom):
			line = content[line_index].split()
			line_index += 1

			### identification
			ai = int(line[0])-1
			molecules[ai] = line[1]
			types[ai] = line[2]

			### assume molecular atom style
			if len(line) == 6 or len(line) == 9:
				r[ai] = line[3:6]

			### assume full atom style
			elif len(line) == 7 or len(line) == 10:
				if i == 0:
					readCharge = True
				charges[ai] = line[3]
				r[ai] = line[4:7]

			### throw error
			else:
				print("Error: Cannot read geometry - unable to surmise atom style.\n")
				sys.exit()

	### get bond information
	bonds = np.zeros((nbond,3), dtype=int)
	if nbond:
		for i in range(len(content)):
			if len(content[i].split()) > 0 and content[i].split()[0] == 'Bonds':
				line_index = i+2
				break
		for i in range(nbond):
			bonds[i,0] = content[line_index].split()[1]
			bonds[i,1] = content[line_index].split()[2]
			bonds[i,2] = content[line_index].split()[3]
			line_index += 1

	### get angle information
	angles = np.zeros((nangle,4), dtype=int)
	if nangle:
		for i in range(len(content)):
			if len(content[i].split()) > 0 and content[i].split()[0] == 'Angles':
				line_index = i+2
				break
		for i in range(nangle):
			angles[i,0] = content[line_index].split()[1]
			angles[i,1] = content[line_index].split()[2]
			angles[i,2] = content[line_index].split()[3]
			angles[i,3] = content[line_index].split()[4]
			line_index += 1


	### get bond information
	if extraLabel is not None:
		line_index = -1
		for i in range(len(content)):
			if len(content[i].split()) > 0 and content[i].split()[0] == extraLabel:
				line_index = i+2
				break

		### check if extras were found
		if line_index == -1:
			print("Error: extra label not found.\n")
			sys.exit()

		nextra = len(content[line_index].split())-1
		extras = np.zeros((natom,nextra))
		for i in range(natom):
			line = content[line_index].split()
			line_index += 1

			ai = int(line[0])-1
			for j in range(nextra):
				extras[ai,j] = line[j+1]

	### results
	output = [ r, molecules, types ]
	if readCharge: output.append(charges)
	output.append(bonds)
	output.append(angles)
	if extraLabel is not None: output.append(extras)
	if getDbox: output.append(dbox3[0])
	if getDbox3: output.append(dbox3)
	return output


### read wham output
def readWham(whamFile, nbin):
	ars.checkFileExist(whamFile, "wham")
	with open(whamFile, 'r') as f:
		content = f.readlines()
	op = np.zeros(nbin)
	PMF = np.zeros(nbin)
	PMF_err = np.zeros(nbin)
	for i in range(nbin):
		op[i] = content[i+2].split()[0]
		PMF[i] = content[i+2].split()[1]
		PMF_err[i] = content[i+2].split()[2]
	return op, PMF, PMF_err


### read my umbrella sampling metadata file
def readUSmeta(metaFile):
	ars.checkFileExist(metaFile, "metadata")
	with open(metaFile, 'r') as f:
		content = f.readlines()
	content = ars.cleanFileContent(content)
	nsim = len(content)-1
	ops = np.zeros(nsim)
	weights = np.zeros(nsim)
	for i in range(nsim):
		ops[i] = float(content[i+1].split()[0])
		weights[i] = float(content[i+1].split()[1])
	return ops, weights


### read cluster file
def readCluster(clusterFile):
	ars.checkFileExist(clusterFile, "cluster")
	with open(clusterFile, 'r') as f:
		content = f.readlines()
	content = ars.cleanFileContent(content)
	ncluster = int(len(content)/2)
	clusters = [None]*ncluster
	for c in range(ncluster):
		indices = content[1+c*2].split()
		clusters[c] = [None]*len(indices)
		for i in range(len(indices)):
			clusters[c][i] = int(indices[i])
	return clusters

	
### read file containing simulation folder names
def readCopies(copiesFile):
	ars.checkFileExist(copiesFile, "copies")
	with open(copiesFile, 'r') as f:
		content = f.readlines()
	content = ars.cleanFileContent(content)
	nsim = len(content)
	copyNames = [None]*nsim
	for i in range(nsim):
		copyNames[i] = content[i].split()[0]
	return copyNames


### remove blank lines and comments from file content
def cleanFileContent(content):
	cleaned = []

	### loop over lines
	for line in content:

		### remove leading and trailing whitespaces
		stripped = line.strip()

		### remove blank lines and whole-line comments
		if not stripped or stripped.startswith("#"):
			continue

		### remove in-line comments
		if "#" in stripped:
			stripped = stripped.split("#", 1)[0].strip()

		### add meaty line to cleaned content
		if stripped:
			cleaned.append(stripped)

	### result
	return cleaned


### search given file for references to arsenal functions
def findArsReferences(searchFile=sys.argv[0], hush=False):

	### read content of search file
	ars.checkFileExist(searchFile, "search")
	with open(searchFile, 'r') as f:
		content = f.read()

	### look for arsenal references in the search file
	pattern = r'\bars\.(\w+)\('
	ars_refs = list(set(re.findall(pattern, content)))
	if "findArsReferences" in ars_refs:
		ars_refs.remove("findArsReferences")

	### print the results
	if not hush:
		print("Arsenal references:")
		for i in range(len(ars_refs)):
			print(f"- {ars_refs[i]}")
	return ars_refs


################################################################################
### File Writers

### write lammps-style geometry
def writeGeo(geoFile, dbox3, r, molecules='auto', types='auto', bonds=None, angles=None, **kwargs):

	### additional keyword args
	natomType	= 'auto'	if 'natomType' not in kwargs else kwargs['natomType']
	nbondType	= 'auto'	if 'nbondType' not in kwargs else kwargs['nbondType']
	nangleType	= 'auto'	if 'nangleType' not in kwargs else kwargs['nangleType']
	masses		= 'auto'	if 'masses' not in kwargs else kwargs['masses']
	charges		= None		if 'charges' not in kwargs else kwargs['charges']
	extras		= None		if 'extras' not in kwargs else kwargs['extras']
	x_precision	= 2			if 'x_precision' not in kwargs else kwargs['x_precision']
	q_precision	= 4			if 'q_precision' not in kwargs else kwargs['q_precision']
	e_precision	= 4			if 'e_precision' not in kwargs else kwargs['e_precision']

	### notes
	# by convention, molecule/type/bond/angle indexing starts at 1, however, there is
	  # nothing in this function that requires this to be the case (it will print the
	  # values it is given).
	# bonds: nbond x 3 (bond type, atom 1, atom2) numpy array
	# angles: nangle x 4 (angle type, atom 1, atom2, atom 3) numpy array

	### count atoms
	natom = len(r)

	### interpret input
	if ars.isnumber(dbox3):
		dbox3 = np.ones(3)*dbox3
	elif not ars.isarray(dbox3) or len(dbox3) != 3:
		print("Flag: Not writing geometry file - dbox3 must be number or 3-element array.")
		return
	if isinstance(molecules, str) and molecules == 'auto':
		molecules = np.zeros(natom, dtype=int)
	if isinstance(types, str) and types == 'auto':
		types = np.ones(natom, dtype=int)
	if bonds is None:
		bonds = np.zeros((0,3), dtype=int)
	if angles is None:
		angles = np.zeros((0,4), dtype=int)

	### inerpret charges
	if charges is None:
		includeCharge = False
	elif isinstance(charges, str) and charges == 'auto':
		includeCharge = True
		charges = np.zeros(natom, dtype=int)
		len_charge = 1
	elif ars.isarray(charges) and len(charges) == natom:
		includeCharge = True
		len_charge = len(str(int(max(charges))))
	else:
		print("Flag: Not writing geometry file - charges must be 'auto' or natom-element array.")
		return

	### inerpret extras
	if extras is None:
		includeExtra = False
	elif ars.isarray(extras) and len(extras) == natom:
		includeExtra = True
		extras = np.array(extras)
		if len(extras.shape) == 1:
			extras = extras.reshape(-1,1)
		nextra = extras.shape[1]
		len_extra = [None]*nextra
		for i in range(nextra):
			len_extra[i] = len(str(int(max(extras[:,i]))))
	else:
		print("Flag: Not writing geometry file - if given, extras must array with natom elements in the first dimension.")
		return

	### count objects
	nmolecule = int(max(molecules))
	nbond = len(bonds)
	nangle = len(angles)

	### interpret input
	if isinstance(natomType, str) and natomType == 'auto':
		natomType = int(max(types))
	elif not ars.isinteger(natomType):
		print("Flag: Not writing geometry file - natomType must be 'auto' or integer.")
		return
	if isinstance(nbondType, str) and nbondType == 'auto':
		if nbond > 0:
			nbondType = int(max(bonds[:,0]))
		else:
			nbondType = 0
	elif not ars.isinteger(nbondType):
		print("Flag: Not writing geometry file - nbondType must be 'auto' or integer.")
		return
	if isinstance(nangleType, str) and nangleType == 'auto':
		if nangle > 0:
			nangleType = int(max(angles[:,0]))
		else:
			nangleType = 0
	elif not ars.isinteger(nangleType):
		print("Flag: Not writing geometry file - nangleType must be 'auto' or integer.")
		return

	### interpret masses
	if isinstance(masses, str) and masses == 'auto':
		masses = np.ones(natomType, dtype=int)
	elif not ars.isarray(masses) or len(masses) != natomType:
		print("Flag: Not writing geometry file - masses must be 'auto' or natomType-element array.")
		return

	### count digits
	len_natom = len(str(natom))
	len_nbond = len(str(nbond))
	len_nangle = len(str(nangle))
	len_nobject = max([len_natom,len_nbond,len_nangle])
	len_nmolecule = len(str(nmolecule))
	len_natomType = len(str(natomType))
	len_nbondType = len(str(nbondType))
	len_nangleType = len(str(nangleType))
	len_nobjectType = max([len_natomType,len_nbondType,len_nangleType])
	len_dbox3 = len(str(int(max(dbox3)/2)))

	### write to file
	with open(geoFile, 'w') as f:

		f.write("## Number of Objects\n")
		f.write(f"\t{natom:<{len_nobject}} atoms\n")
		if nbond:
			f.write(f"\t{nbond:<{len_nobject}} bonds\n")
		if nangle:
			f.write(f"\t{nangle:<{len_nobject}} angles\n")

		f.write("\n## Number of Object Types\n")
		f.write(f"\t{natomType:<{len_nobjectType}} atom types\n")
		if nbondType:
			f.write(f"\t{nbondType:<{len_nobjectType}} bond types\n")
		if nangleType:
			f.write(f"\t{nangleType:<{len_nobjectType}} angle types\n")

		f.write("\n## Simulation Box\n")
		f.write(f"\t{-dbox3[0]/2:>{len_dbox3+x_precision+2}.{x_precision}f} {dbox3[0]/2:>{len_dbox3+x_precision+1}.{x_precision}f} xlo xhi\n")
		f.write(f"\t{-dbox3[1]/2:>{len_dbox3+x_precision+2}.{x_precision}f} {dbox3[1]/2:>{len_dbox3+x_precision+1}.{x_precision}f} ylo yhi\n")
		f.write(f"\t{-dbox3[2]/2:>{len_dbox3+x_precision+2}.{x_precision}f} {dbox3[2]/2:>{len_dbox3+x_precision+1}.{x_precision}f} zlo zhi\n")

		f.write("\nMasses\n\n")
		for i in range(natomType):
			f.write(f"\t{i+1:<{len_natomType}} {masses[i]}\n")

		f.write("\nAtoms\n\n")
		for i in range(natom):
			f.write(f"\t{i+1:<{len_natom}}" + \
					f" {int(molecules[i]):<{len_nmolecule}}" + \
					f" {int(types[i]):<{len_natomType}}")
			if includeCharge:
				f.write(f" {charges[i]:>{len_charge+q_precision+1}.{q_precision}f}") 
			f.write(f"  {r[i,0]:>{len_dbox3+x_precision+2}.{x_precision}f}" + \
					 f" {r[i,1]:>{len_dbox3+x_precision+2}.{x_precision}f}" + \
					 f" {r[i,2]:>{len_dbox3+x_precision+2}.{x_precision}f}\n")

		if nbond:
			f.write("\nBonds\n\n")
			for i in range(nbond):
				f.write(f"\t{i+1:<{len_nbond}} " + \
						f"{int(bonds[i,0]):<{len_nbondType}}  " + \
						f"{int(bonds[i,1]):<{len_natom}} " + \
						f"{int(bonds[i,2]):<{len_natom}}\n")

		if nangle:
			f.write("\nAngles\n\n")
			for i in range(nangle):
				f.write(f"\t{i+1:<{len_nangle}} " + \
						f"{int(angles[i,0]):<{len_nangleType}}  " + \
						f"{int(angles[i,1]):<{len_natom}} " + \
						f"{int(angles[i,2]):<{len_natom}} " + \
						f"{int(angles[i,3]):<{len_natom}}\n")

		if includeExtra:
			f.write("\nExtras\n\n")
			for i in range(natom):
				f.write(f"\t{i+1:<{len_natom}}")
				for j in range(nextra):
					f.write(f" {extras[i][j]:>{len_extra[j]+e_precision+2}.{e_precision}f}")
				f.write("\n")


### load array of variables from pickle file
def unpickle(pklFile, dims=None, hushExtraFlag=False):

	### read file
	ars.checkFileExist(pklFile, "pickle")
	with open(pklFile, 'rb') as f:
		cucumber = pickle.load(f)

	### check and trim content
	if dims is not None:
		nvar = len(dims)

		### make into array if single value
		if not ars.isarray(cucumber):
			cucumber = [ cucumber ]

		### make sure pickle array is long enough
		if len(cucumber) < len(dims):
			print("Error: pickle array contained fewer elements than expected.\n")
			sys.exit()

		### trim pickle array if too long
		elif len(cucumber) > len(sims):
			if not hushExtraFlag:
				print("Flag: pickle array contained more elements than expected, ignoring some elements of the pickle array.")
			cucumber = cucumber[:len(dims)]

		### check variables
		for i in range(nvar):

			### unknown dimension
			if dims[i] == None:
				continue

			### single value expected
			elif dims[i] == 0:
				if not ars.isnumber(cucumber[i]):
					print(f"Error: element {i} of the pickle array does not match the expected data type (number).\n")
					sys.exit()

			### if not single value, must be array
			elif not ars.isarray(cucumber[i]):
				print(f"Error: element {i} of the pickle array does not match the expected data type (array).\n")
				sys.exit()

			### if numpy array, must have correct shape
			elif isinstance(cucumber[i], np.ndarray):
				if len(cucumber[i].shape) != dims[i]:
					print(f"Error: element {i} of the pickle array does not have the expected shape ({dims[i]}).\n")
					sys.exit()

	### result
	return cucumber


### create arsenal file that contains only the functions necessary for the scripts in a given folder
def deployArsenal(srcFold=os.getcwd()+"/"):

	### location of arsenal file and deployed arsenal file
	arsFile = "/Users/dduke/Programs/arsenal/arsenal.py"
	arsDepFile = srcFold + "armament.py"

	### get arsenal code
	checkFileExist(arsFile, 'arsenal')
	with open(arsFile, 'r') as f:
		ars_code = f.readlines()

	### gather all arsenal function definitions
	ars_tree = ast.parse("".join(ars_code))
	ars_functions = {}
	for node in ast.walk(ars_tree):
		if isinstance(node, ast.FunctionDef) and node.name != "deployArsenal":
			start_line = node.lineno - 2
			end_line = max(child.end_lineno for child in ast.walk(node) if hasattr(child, 'end_lineno'))
			ars_functions[node.name] = (start_line, end_line)

	### find all arsenal references
	ars_refs = []    
	for root, dirs, files in os.walk(srcFold):
		dirs[:] = [d for d in dirs if not d.startswith(".") and d != "__pycache__"]
		for file in files:
			if file.endswith(".py"):
				srcFile = os.path.join(root, file)
				ars_refs.extend(ars.findArsReferences(srcFile, hush=True))

	### write deployed arsenal file
	deployed_functions = []
	with open(arsDepFile, 'w') as f:
		for func in ars_functions.keys():
			if func in ars_refs:
				deployed_functions.append(func)
				start, end = ars_functions[func]
				f.writelines(ars_code[start:end])
				f.write("\n\n")

	### check for self-references
	while True:
		ars_refs = ars.findArsReferences(arsDepFile, hush=True)
		ars_refs = [item for item in ars_refs if item not in set(deployed_functions)]
		if len(ars_refs) == 0:
			break
		with open(arsDepFile, 'a') as f:
			for func in ars_functions.keys():
				if func in ars_refs:
					deployed_functions.append(func)
					start, end = ars_functions[func]
					f.writelines(ars_code[start:end])
					f.write("\n\n")

	### read contents of deployed arsenal file
	with open(arsDepFile, 'r') as f:
		ars_dep_code = f.readlines()
	ars_dep_tree = ast.parse("".join(ars_dep_code))

	### extract all function references in deployed arsenal
	all_refs = set()
	for node in ast.walk(ars_dep_tree):
		if isinstance(node, ast.Name):
			all_refs.add(node.id)

	### extract initial imports in arsenal code
	imports = []
	for i, line in enumerate(ars_code):
		if line.strip() == "":
			break
		if line.strip().startswith("import") or line.strip().startswith("from"):
			if line.strip() == "import arsenal as ars":
				imports.append("import armament as ars")
			else:
				imports.append(line.strip())

	### compile list of imports referenced in deployed arsenal
	used_imports = []
	for imp in imports:

		### handle "from module import ..."
		if imp.startswith("from"):
			module, _, items = imp.partition("import")
			module = module.strip()[5:]  # Remove "from "
			items = [item.strip() for item in items.split(",")]
			used_items = [item for item in items if item.split()[0] in all_refs]
			if used_items:
				used_imports.append(f"from {module} import {', '.join(used_items)}")

		### handle "import module" or "import module as alias"
		elif imp.startswith("import"):
			parts = imp[7:].strip().split(" as ")
			module = parts[0].strip()
			alias = parts[1].strip() if len(parts) > 1 else module
			if module in all_refs or alias in all_refs or any(alias + "." in name or module + "." in name for name in all_refs):
				used_imports.append(imp)

	### re-write deployed arsenal with the used imports
	with open(arsDepFile, 'w') as f:
		f.writelines(line + "\n" for line in used_imports)
		f.write("\n")
		f.writelines(ars_dep_code)


################################################################################
### Plotting

### set pretty matplotlib defaults
def magicPlot(pubReady=False):

	### set default magic settings
	params = {
		'figure.figsize'	: '8,6',
		'font.family'		: 'Times',
		'text.usetex'		: True,
		'errorbar.capsize'	: 3,
		'lines.markersize'	: 6,
		'legend.fontsize'	: 14,
		'xtick.labelsize'	: 14,
		'ytick.labelsize'	: 14,
		'axes.labelsize'	: 14,
		'axes.titlesize'	: 16,
		'axes.axisbelow'	: True,
		'grid.color'		: (0.5, 0.5, 0.5, 0.3)
	}
	plt.rcParams.update(params)

	### set default magic settings
	paramsPub = {
		'figure.figsize'	: '8,6',
		'font.family'		: 'Times',
		'text.usetex'		: True,
		'errorbar.capsize'	: 3,
		'lines.markersize'	: 6,
		'legend.fontsize'	: 20,
		'xtick.labelsize'	: 20,
		'ytick.labelsize'	: 20,
		'axes.labelsize'	: 20,
		'axes.titlesize'	: 24,
		'axes.axisbelow'	: True,
		'grid.color'		: (0.5, 0.5, 0.5, 0.3)
	}

	if not pubReady:
		plt.rcParams.update(params)
	else:
		plt.rcParams.update(paramsPub)



### plot the convergence of a varaible
def plotConv(A, Alabel=None, title=None, figLabel='auto', Alim='auto', Xlim='auto', Xlabel='auto', dt_per_use='auto', **kwargs):

	### additional keyword args
	ax = None if 'ax' not in kwargs else kwargs['ax']

	### interpret input
	if isinstance(Alim, str) and Alim == 'auto':
		Amin = min(A)
		Amax = max(A)
		Alim = [ Amin-(Amax-Amin)*0.6, Amax+(Amax-Amin)*0.4 ]
	elif not ars.isarray(Alim) or len(Alim) != 2:
		print("Flag: Skipping convergence plot - varaible limits must be either 'auto' or 2-element array.")
		return
	if isinstance(Xlim, str) and Xlim == 'auto':
		Xlim = [ 0, len(A)-1 ]
	elif not ars.isarray(Xlim) or len(Xlim) != 2:
		print("Flag: Skipping convergence plot - x-axis limits must be either 'auto' or 2-element array.")
		return
	if isinstance(dt_per_use, str) and dt_per_use == 'auto':
		time = np.arange(len(A))
	elif ars.isnumber(dt_per_use):
		time = np.arange(len(A))*dt_per_use
	else:
		print("Flag: Skipping convergence plot - time step must be either 'auto' or number.")
		return
	if ax is not None:
		if figLabel != 'auto' and figLabel is not None:
			print("Warining: Unused figure label input.")
			figLabel = None
	elif figLabel == 'auto':
		figLabel = "Conv"

	### calculate mean and sem
	avg = np.zeros(len(A))
	sem = np.zeros(len(A))
	for i in range(len(A)):
		avg[i] = np.mean(A[0:i+1])
		if i > 0:
			sem[i] = ars.calcSEMautocorr(A[0:i+1], hush=True)

	### initialize figure
	if ax is not None:
		plt.sca(ax)
	else:
		plt.figure(figLabel)

	### plot data, mean, sem
	plt.scatter(time, A, color='green')
	plt.plot(time, avg, color='purple')
	plt.fill_between(time, avg+sem, avg-sem, color='purple', alpha=0.3, linewidth=0)
	plt.xlim(Xlim)
	plt.ylim(Alim)
	if Xlabel == 'auto':
		if dt_per_use == 'auto':
			plt.xlabel("Time Step")
		else:
			plt.xlabel("Time")
	elif Xlabel is not None:
		plt.xlabel(Xlabel)
	if Alabel is not None:
		plt.ylabel(Alabel)
	plt.legend(['Data','Mean','SEM'], loc='lower right')
	if title is not None:
		plt.title(title)


### plot a some nice points
def plotPoints(X, Y, title=None, figLabel="Points", Xlim='auto', Ylim='auto', Xlabel=None, Ylabel=None, **kwargs):
	
	### additional keyword args
	S			= None		if 'S' not in kwargs else kwargs['S']
	marker		= 'o'		if 'marker' not in kwargs else kwargs['marker']
	color		= 'black'	if 'color' not in kwargs else kwargs['color']
	edgecolor	= None		if 'edgecolor' not in kwargs else kwargs['edgecolor']
	alpha		= None		if 'alpha' not in kwargs else kwargs['alpha']
	zorder		= None		if 'zorder' not in kwargs else kwargs['zorder']
	label		= None		if 'label' not in kwargs else kwargs['label']
	E			= None		if 'E' not in kwargs else kwargs['E']
	ecolor		= None		if 'ecolor' not in kwargs else kwargs['ecolor']
	ecapsize	= None		if 'ecapsize' not in kwargs else kwargs['ecapsize']
	elinewidth	= None		if 'elinewidth' not in kwargs else kwargs['elinewidth']
	ax			= None		if 'ax' not in kwargs else kwargs['ax']

	### interpret input
	if ax is not None:
		if figLabel != 'auto' and figLabel is not None:
			print("Warining: Unused figure label input.")
			figLabel = None
	elif figLabel == 'auto':
		figLabel = "Points"

	### initialize figure
	if ax is not None:
		plt.sca(ax)
	else:
		plt.figure(figLabel)

	### adjust markersize for scatter
	if S is not None:
		S = S**2 

	### plot errorbars
	if E is not None:
		plt.errorbar(X, Y, E, fmt='none', ecolor=ecolor, capsize=ecapsize, linewidth=elinewidth, capthick=elinewidth)

	### plot points
	plt.scatter(X, Y, S, marker=marker, color=color, edgecolor=edgecolor, alpha=alpha, zorder=zorder, label=label)

	### format
	if Xlim != 'auto':
		plt.xlim(Xlim)
	if Ylim != 'auto':
		plt.ylim(Ylim)
	if Xlabel is not None:
		plt.xlabel(Xlabel)
	if Ylabel is not None:
		plt.ylabel(Ylabel)
	if title is not None:
		plt.title(title)


### plot a nice line
def plotLine(X, Y, title=None, figLabel='auto', Xlim='auto', Ylim='auto', Xlabel=None, Ylabel=None, **kwargs):

	### additional keyword args
	color		= None		if 'color' not in kwargs else kwargs['color']
	linestyle	= '-'		if 'linestyle' not in kwargs else kwargs['linestyle']
	marker		= None		if 'marker' not in kwargs else kwargs['marker']
	markersize	= None		if 'markersize' not in kwargs else kwargs['markersize']
	alpha		= None		if 'alpha' not in kwargs else kwargs['alpha']
	zorder		= None		if 'zorder' not in kwargs else kwargs['zorder']
	label		= None		if 'label' not in kwargs else kwargs['label']
	E			= None		if 'E' not in kwargs else kwargs['E']
	ecolor		= None		if 'ecolor' not in kwargs else kwargs['ecolor']
	ecapsize	= None		if 'ecapsize' not in kwargs else kwargs['ecapsize']
	elinewidth	= None		if 'elinewidth' not in kwargs else kwargs['elinewidth']
	label		= None		if 'label' not in kwargs else kwargs['label']
	ax			= None		if 'ax' not in kwargs else kwargs['ax']

	### interpret input
	if ax is not None:
		if figLabel != 'auto' and figLabel is not None:
			print("Warining: Unused figure label input.")
			figLabel = None
	elif figLabel == 'auto':
		figLabel = "Line"

	### initialize figure
	if ax is not None:
		plt.sca(ax)
	else:
		plt.figure(figLabel)
	
	### plot errorbars
	if E is not None:
		plt.errorbar(X, Y, E, fmt='none', ecolor=ecolor, capsize=ecapsize, linewidth=elinewidth, capthick=elinewidth)

	### plot line
	plt.plot(X, Y, color=color, linestyle=linestyle, marker=marker, markersize=markersize, alpha=alpha, zorder=zorder,label=label)

	### format
	if Xlim != 'auto':
		plt.xlim(Xlim)
	if Ylim != 'auto':
		plt.xlim(Ylim)
	if Xlabel is not None:
		plt.xlabel(Xlabel)
	if Ylabel is not None:
		plt.ylabel(Ylabel)
	if title is not None:
		plt.title(title)


### plot a nice histogram
def plotHist(A, Alabel=None, title=None, figLabel='auto', nbin='auto', Alim_bin='auto', Alim_plot='auto', Ylabel='auto', useDensity=False, rng=np.random, **kwargs):

	### additional keyword args
	weights			= None		if 'weights' not in kwargs else kwargs['weights']
	assumeIID		= False		if 'assumeIID' not in kwargs else kwargs['assumeIID']
	plotBins		= True		if 'plotBins' not in kwargs else kwargs['plotBins']
	plotSteps		= False		if 'plotSteps' not in kwargs else kwargs['plotSteps']
	plotLine		= False		if 'plotLine' not in kwargs else kwargs['plotLine']
	plotGauss		= False		if 'plotGauss' not in kwargs else kwargs['plotGauss']
	plotNorm		= False		if 'plotNorm' not in kwargs else kwargs['plotNorm']
	plotLogNorm		= False		if 'plotLogNorm' not in kwargs else kwargs['plotLogNorm']
	color			= None		if 'color' not in kwargs else kwargs['color']
	edgecolor		= 'auto'	if 'edgecolor' not in kwargs else kwargs['edgecolor']
	alpha			= 0.6		if 'alpha' not in kwargs else kwargs['alpha']
	solidify		= False		if 'solidify' not in kwargs else kwargs['solidify']
	label			= None		if 'label' not in kwargs else kwargs['label']
	gauss_label		= 'auto'	if 'gauss_label' not in kwargs else kwargs['gauss_label']
	gauss_color		= 'match'	if 'gauss_color' not in kwargs else kwargs['gauss_color']
	norm_label		= 'auto'	if 'norm_label' not in kwargs else kwargs['norm_label']
	norm_color		= 'match'	if 'norm_color' not in kwargs else kwargs['norm_color']
	logNorm_label	= 'auto'	if 'logNorm_label' not in kwargs else kwargs['logNorm_label']
	logNorm_color	= 'match'	if 'logNorm_color' not in kwargs else kwargs['logNorm_color']
	plotErrSR		= False		if 'plotErrSR' not in kwargs else kwargs['plotErrSR']
	plotErrBS		= False		if 'plotErrBS' not in kwargs else kwargs['plotErrBS']
	plotGaussErr	= False		if 'plotGaussErr' not in kwargs else kwargs['plotGaussErr']
	plotNormErr		= False		if 'plotNormErr' not in kwargs else kwargs['plotNormErr']
	plotLogNormErr	= False		if 'plotLogNormErr' not in kwargs else kwargs['plotLogNormErr']
	sigma_ci		= 1			if 'sigma_ci' not in kwargs else kwargs['sigma_ci']
	ecolor			= 'auto'	if 'ecolor' not in kwargs else kwargs['ecolor']
	plotAvg			= False		if 'plotAvg' not in kwargs else kwargs['plotAvg']
	plotMed			= False		if 'plotMed' not in kwargs else kwargs['plotMed']
	plotStd			= False		if 'plotStd' not in kwargs else kwargs['plotStd']
	avg_label		= 'auto'	if 'avg_label' not in kwargs else kwargs['avg_label']
	avg_color		= 'red'		if 'avg_color' not in kwargs else kwargs['avg_color']
	med_label		= 'auto'	if 'med_label' not in kwargs else kwargs['med_label']
	med_color		= 'red'		if 'med_color' not in kwargs else kwargs['med_color']
	std_label		= 'auto'	if 'std_label' not in kwargs else kwargs['std_label']
	std_color		= 'red'		if 'std_color' not in kwargs else kwargs['std_color']
	ax				= None		if 'ax' not in kwargs else kwargs['ax']
	hush			= False		if 'hush' not in kwargs else kwargs['hush']

	### determine whether plotting raw data
	plotData = plotBins or plotSteps or plotLine

	### ensure valid input
	if not (plotBins or plotSteps or plotLine or plotGauss or plotNorm or plotLogNorm):
		print("Flag: Skipping histogram plot - at least one plot type (bins, steps, line, gauss, normal, log normal) must be true.")
		return
	if weights is not None and any(weights<0):
		print("Warning: Skipping histogram plot - weights must be non-negative.")
		return
	if sigma_ci <= 0:
		print("Flag: Skipping histogram plot - confidence interval sigma must be positive.")
		return
	if (plotErrSR or plotErrBS) and not plotData:
		print("Flag: Removing histogram error bars - at least one data plot (bins, steps, line) must be true to plot error bars.")
		plotErrSR, plotErrBS = False, False
	if plotLogNorm and any(A<=0):
		print("Flag: Removing histogram log normal fit - data must be positive.")
		plotLogNorm = False
	if plotGaussErr and not plotGauss:
		print("Flag: Removing histogram Gaussian KDE error bars - cannot plot error bars without plotting fit.")
		plotGaussErr = False
	if plotNormErr and not plotNorm:
		print("Flag: Removing histogram normal fit error bars - cannot plot error bars without plotting fit.")
		plotNormErr = False
	if plotLogNormErr and not plotLogNorm:
		print("Flag: Removing histogram log normal error bars - cannot plot error bars without plotting fit.")
		plotLogNormErr = False

	### interpret input
	if isinstance(nbin, str) and nbin == 'auto':
		nbin = ars.optbins(A, 50)
	elif not ars.isinteger(nbin):
		print("Flag: Skipping histogram plot - number of histogram bins must be either 'auto' or integer.")
		return
	if isinstance(Alim_bin, str) and Alim_bin == 'auto':
		Alim_bin = [ min(A), max(A) ]
		if Alim_bin[0] == Alim_bin[1]:
			print("Flag: Skipping hitogram plot - all values are the same.")
			return
	elif not ars.isarray(Alim_bin) or len(Alim_bin) != 2:
		print("Flag: Skipping histogram plot - variable limits must be either 'auto' or 2-element array.")
		return
	if isinstance(Alim_plot, str) and Alim_plot == 'auto':
		dAbin = (Alim_bin[1]-Alim_bin[0])/nbin
		Alim_plot = [ Alim_bin[0]-dAbin/2, Alim_bin[1]+dAbin/2 ]
	elif not ars.isarray(Alim_plot) or len(Alim_plot) != 2:
		if Alim_plot is not None:
			print("Flag: Skipping histogram plot - variable limits must be None, 'auto', or 2-element array.")
			return
	if isinstance(edgecolor, str) and edgecolor == 'auto':
		edgecolor = 'black'
		if plotSteps or plotLine or plotGauss or plotNorm or plotLogNorm:
			edgecolor = None
	if isinstance(ecolor, str) and ecolor == 'auto':
		ecolor = edgecolor if edgecolor is not None else color
	if ax is not None:
		if figLabel != 'auto' and figLabel is not None:
			print("Flag: Unused figure label input.")
			figLabel = None
	elif figLabel == 'auto':
		figLabel = "Hist"

	### check for unused labels
	if not plotGauss:
		if not (gauss_label is None or gauss_label == 'auto'):
			print("Flag: Unused label input for Gaussian KDE.")
		gauss_label = None
	if not plotNorm:
		if not (norm_label is None or norm_label == 'auto'):
			print("Flag: Unused label input for normal fit.")
		norm_label = None
	if not plotLogNorm:
		if not (logNorm_label is None or logNorm_label == 'auto'):
			print("Flag: Unused label input for log normal fit.")
		logNorm_label = None
	if not plotAvg:
		if not (avg_label is None or avg_label == 'auto'):
			print("Flag: Unused label input for average line.")
		avg_label = None
	if not plotMed:
		if not (med_label is None or med_label == 'auto'):
			print("Flag: Unused label input for median line.")
		med_label = None
	if not plotStd:
		if not (std_label is None or std_label == 'auto'):
			print("Flag: Unused label input for standard deviation lines.")
		std_label = None

	### interpret fit labels
	if gauss_label == 'auto':
		if label is not None:
			gauss_label = None if plotData else 'match'
		else:
			gauss_label = "Gaussian KDE" if plotData else None
	if norm_label == 'auto':
		if label is not None:
			norm_label = None if plotData else 'match'
		else:
			norm_label = "Normal Fit" if plotData else None
	if logNorm_label == 'auto':
		if label is not None:
			logNorm_label = None if plotData else 'match'
		else:
			logNorm_label = "Log Normal Fit" if plotData else None

	### interpret line labels
	if avg_label == 'auto':
		avg_label = f"$\\mu$"
	elif avg_label is not None:
		avg_label = "$\\mu_{\\textrm{" + f"{avg_label}" + "}}$"
	if med_label == 'auto':
		med_label = f"$M$"
	elif med_label is not None:
		med_label = "$M_{\\textrm{" + f"{med_label}" + "}}$"
	if std_label == 'auto':
		std_label = f"$\\sigma$"
	elif std_label is not None:
		std_label = "$\\sigma_{\\textrm{" + f"{std_label}" + "}}$"

	### determine whether plotting legend
	legend = False
	if label or gauss_label or norm_label or logNorm_label or avg_label or med_label or std_label:
		legend = True

	### numpify data
	A = np.asarray(A)
	if weights is not None:
		weights = np.asarray(weights)

	### calculate histogram
	heights, edges = np.histogram(A, nbin, weights=weights, range=Alim_bin)
	centers = edges[:len(edges)-1] + 1/2*(edges[1]-edges[0])
	width_bin = (Alim_bin[1]-Alim_bin[0])/nbin

	### determine how much of the data is plotted
	area_plot = sum(heights)*width_bin
	if weights is not None:
		area_full = sum(weights)*width_bin
	else:
		area_full = len(A)*width_bin

	### determine density scaling
	scale = area_plot if useDensity else 1

	### calculate statistics
	if plotAvg or plotStd:
		avg = np.mean(A)
	if plotMed:
		med = np.median(A)
	if plotStd:
		std = np.std(A)

	### assume weighted data is uncorrelated
	if weights is not None:
		assumeIID = True

	### calculate number of independent data points
	if plotErrSR or plotErrBS or plotGaussErr or plotNormErr or plotLogNormErr:
		if assumeIID:
			tau_int = 1
		else:
			tau_int = 1 + 2*ars.calcCorrTime(A)
			if not hush and tau_int > 1:
				print(f"Using integrated correlation time estimate {tau_int:0.2f} for error bars.")

	### fit parameters
	npoint_fit = 100
	nbootstrap = 1000
	percentile_lower = norm.cdf(-sigma_ci)*100
	percentile_upper = norm.cdf(sigma_ci)*100

	### initialize figure
	if ax is not None:
		plt.sca(ax)
	else:
		plt.figure(figLabel)

	### initialize label flag
	isDataLabeled = True if label is None else False

	### plot solid bin background
	if (plotBins or plotSteps) and solidify:
		plt.bar(centers, heights/scale, width_bin, color='white')

	### plot data as bins
	if plotBins:		
		bars = plt.bar(centers, heights/scale, width_bin, color=color, alpha=alpha, edgecolor=edgecolor)
		color = bars[0].get_facecolor()
		if not isDataLabeled:
			bars[0].set_label(label)
			isDataLabeled = True

	### plot data as steps
	if plotSteps:
		bars = plt.hist(A, nbin, weights=weights, range=Alim_bin, color=color, linewidth=2, histtype='step')[2]
		color = bars[0].get_edgecolor()
		if not isDataLabeled:
			bars[0].set_label(label)
			isDataLabeled = True

	### plot data as line
	if plotLine:
		line = plt.plot(centers, heights, color=color, linewidth=2)[0]
		color = line.get_color()
		if not isDataLabeled:
			line.set_label(label)
			isDataLabeled = True

	### square root error bars on data
	if plotErrSR:
		Y_err = np.sqrt(tau_int*heights)
		plt.errorbar(centers, heights/scale, Y_err/scale, fmt='none', linewidth=1, color=ecolor, alpha=alpha)

	### bootstrap error bars on data
	if plotErrBS:
		boot = np.zeros((nbootstrap, nbin))
		probs = weights / np.sum(weights) if weights is not None else None
		for b in range(nbootstrap):
			idxs = rng.choice(len(A), size=len(A), replace=True, p=probs)
			boot[b] = np.histogram(A[idxs], nbin, range=Alim_bin)[0]
		heights_lower = np.percentile(boot, percentile_lower, axis=0)
		heights_upper = np.percentile(boot, percentile_upper, axis=0)
		heights_err = np.vstack((heights-heights_lower, heights_upper-heights))
		heights_err *= np.sqrt(tau_int)
		plt.errorbar(centers, heights/scale, heights_err/scale, fmt='none', linewidth=1, color=ecolor, alpha=alpha)

	### gaussian kernel-density estimate
	if plotGauss:
		X = np.linspace(Alim_bin[0], Alim_bin[1], npoint_fit)
		Y = gaussian_kde(A, weights=weights)(X)*area_full
		curve = plt.plot(X, Y/scale, color=color, linewidth=2, alpha=1)[0]
		color = curve.get_color()
		if gauss_color != 'match':
			curve.set_color(gauss_color)
		if gauss_label == 'match':
			if not isDataLabeled:
				curve.set_label(label)
				isDataLabeled = True
		elif gauss_label is not None:
			curve.set_label(gauss_label)

		### gaussian KDE error
		if plotGaussErr:
			boot = np.zeros((nbootstrap, npoint_fit))
			probs = weights / np.sum(weights) if weights is not None else None
			for b in range(nbootstrap):
				idxs = rng.choice(len(A), size=len(A), replace=True, p=probs)
				boot[b] = gaussian_kde(A[idxs])(X)*area_full
			Y_lower = np.percentile(boot, percentile_lower, axis=0)
			Y_upper = np.percentile(boot, percentile_upper, axis=0)
			Y_lower = Y - (Y-Y_lower)*np.sqrt(tau_int)
			Y_upper = Y + (Y_upper-Y)*np.sqrt(tau_int)
			patch = plt.fill_between(X, Y_lower/scale, Y_upper/scale, color=color, alpha=0.2, linewidth=0)
			if gauss_color != 'match':
				patch.set_color(gauss_color)

	### normal fit
	if plotNorm:
		mu, sigma = ars.calcNormStats(A, weights)
		X = np.linspace(Alim_bin[0], Alim_bin[1], npoint_fit)
		Y = norm.pdf(X, loc=mu, scale=sigma)*area_full
		curve = plt.plot(X, Y/scale, color=color, linewidth=2, alpha=1)[0]
		color = curve.get_color()
		if norm_color != 'match':
			curve.set_color(norm_color)
		if norm_label == 'match':
			if not isDataLabeled:
				curve.set_label(label)
				isDataLabeled = True
		elif norm_label is not None:
			curve.set_label(norm_label)

		### normal fit error
		if plotNormErr:
			probs = weights / np.sum(weights) if weights is not None else None
			boot = np.zeros((nbootstrap, npoint_fit))
			for b in range(nbootstrap):
				idxs = rng.choice(len(A), size=Neff, replace=True, p=probs)
				mu, sigma = ars.calcNormStats(A[idxs])
				boot[b] = norm.pdf(X, loc=mu, scale=sigma)*area_full
			Y_lower = np.percentile(boot, percentile_lower, axis=0)
			Y_upper = np.percentile(boot, percentile_upper, axis=0)
			Y_lower = Y - (Y-Y_lower)*np.sqrt(tau_int)
			Y_upper = Y + (Y_upper-Y)*np.sqrt(tau_int)
			patch = plt.fill_between(X, Y_lower/scale, Y_upper/scale, color=color, alpha=0.2, linewidth=0)
			if norm_color != 'match':
				patch.set_color(norm_color)

	### log normal fit
	if plotLogNorm:
		logA = np.log(A)
		mu, sigma = ars.calcNormStats(logA, weights)
		X = np.linspace(Alim_bin[0], Alim_bin[1], npoint_fit)
		Y = lognorm.pdf(X, s=sigma, scale=np.exp(mu))*area_full
		curve = plt.plot(X, Y/scale, color=color, linewidth=2, alpha=1)[0]
		color = curve.get_color()
		if logNorm_color != 'match':
			curve.set_color(logNorm_color)
		if logNorm_label == 'match':
			if not isDataLabeled:
				curve.set_label(label)
				isDataLabeled = True
		elif logNorm_label is not None:
			curve.set_label(logNorm_label)

		### log normal fit error
		if plotLogNormErr:
			probs = weights / np.sum(weights) if weights is not None else None
			boot = np.zeros((nbootstrap, npoint_fit))
			for b in range(nbootstrap):
				idxs = rng.choice(len(A), size=Neff, replace=True, p=probs)
				mu, sigma = ars.calcNormStats(logA[idxs])
				boot[b] = lognorm.pdf(X, s=sigma, scale=np.exp(mu))*area_full
			Y_lower = np.percentile(boot, percentile_lower, axis=0)
			Y_upper = np.percentile(boot, percentile_upper, axis=0)
			Y_lower = Y - (Y-Y_lower)*np.sqrt(tau_int)
			Y_upper = Y + (Y_upper-Y)*np.sqrt(tau_int)
			patch = plt.fill_between(X, Y_lower/scale, Y_upper/scale, color=color, alpha=0.2, linewidth=0)
			if logNorm_color != 'match':
				patch.set_color(logNorm_color)

	### average line
	if plotAvg:
		line = plt.axvline(avg, linestyle='--', linewidth=2)
		if avg_color == 'match':
			line.set_color(color)
		else:
			line.set_color(avg_color)
		if avg_label is not None:
			line.set_label(f"{avg_label} = {avg:0.2f}")

	### median line
	if plotMed:
		line = plt.axvline(med, linestyle='-.', linewidth=2)
		if med_color == 'match':
			line.set_color(color)
		else:
			line.set_color(med_color)
		if med_label is not None:
			line.set_label(f"{med_label} = {med:0.2f}")

	### standard deviation lines
	if plotStd:
		line = plt.axvline(avg+std, linestyle=':', linewidth=2)
		if std_color == 'match':
			line.set_color(color)
		else:
			line.set_color(std_color)
		line = plt.axvline(avg-std, linestyle=':', linewidth=2, alpha=1)
		if std_color == 'match':
			line.set_color(color)
		else:
			line.set_color(std_color)
		if std_label is not None:
			line.set_label(f"{std_label} = {std:0.2f}")

	### formatting
	plt.ylim(bottom=0, auto=True)
	if Alim_plot is not None:
		plt.xlim(Alim_plot)
	if Alabel is not None:
		plt.xlabel(Alabel)
	if Ylabel == 'auto':
		if useDensity:
			plt.ylabel("Density")
		else:
			plt.ylabel("Count")
	elif Ylabel is not None:
		plt.ylabel(Ylabel)
	if title is not None:
		plt.title(title)
	if legend:
		plt.legend()


### plot a nice 2D histogram
def plotHist2D(A, B, Alabel=None, Blabel=None, title=None, figLabel='auto', nbin='auto', Alim_bin='auto', Blim_bin='auto', Alim_plot='auto', Blim_plot='auto', useDensity=False, **kwargs):

	### additional keyword args
	ax = None if 'ax' not in kwargs else kwargs['ax']

	### interpret input
	if isinstance(nbin, str) and nbin == 'auto':
		nbin = ars.optbins(A, 50)
	elif not ars.isinteger(nbin):
		print("Flag: Skipping histogram plot - number of histogram bins must be either 'auto' or integer.")
		return
	if isinstance(Alim_bin, str) and Alim_bin == 'auto':
		Alim_bin = [ min(A), max(A) ]
	elif not ars.isarray(Alim_bin) or len(Alim_bin) != 2:
		print("Flag: Skipping histogram plot - variable limits must be either 'auto' or 2-element array.")
		return
	if isinstance(Blim_bin, str) and Blim_bin == 'auto':
		Blim_bin = [ min(B), max(B) ]
	elif not ars.isarray(Blim_bin) or len(Blim_bin) != 2:
		print("Flag: Skipping histogram plot - variable limits must be either 'auto' or 2-element array.")
		return
	if isinstance(Alim_plot, str) and Alim_plot == 'auto':
		dAbin = (Alim_bin[1]-Alim_bin[0])/nbin
		Alim_plot = [ Alim_bin[0], Alim_bin[1] ]
	elif not ars.isarray(Alim_plot) or len(Alim_plot) != 2:
		if Alim_plot is not None:
			print("Flag: Skipping histogram plot - variable limits must be None, 'auto', or 2-element array.")
			return
	if isinstance(Blim_plot, str) and Blim_plot == 'auto':
		dBbin = (Blim_bin[1]-Blim_bin[0])/nbin
		Blim_plot = [ Blim_bin[0], Blim_bin[1] ]
	elif not ars.isarray(Blim_plot) or len(Blim_plot) != 2:
		if Blim_plot is not None:
			print("Flag: Skipping histogram plot - variable limits must be None, 'auto', or 2-element array.")
			return
	if ax is not None:
		if figLabel != 'auto' and figLabel is not None:
			print("Warining: Unused figure label input.")
			figLabel = None
	elif figLabel == 'auto':
		figLabel = "Hist2D"

	### initialize figure
	if ax is not None:
		plt.sca(ax)
	else:
		plt.figure(figLabel)

	### plot histogram
	plt.hist2d(A, B, nbin, range=[Alim_bin,Blim_bin], density=useDensity)
	if Alim_plot is not None:
		plt.xlim(Alim_plot)
	if Blim_plot is not None:
		plt.ylim(Blim_plot)
	if Alabel is not None:
		plt.xlabel(Alabel)
	if Blabel is not None:
		plt.ylabel(Blabel)
	if title is not None:
		plt.title(title)


### calculate and plot PMF from histogram
def plotPMF(A, Alabel=None, title=None, figLabel='auto', nbin='auto', Alim_bin='auto', Alim_plot='auto', **kwargs):

	### additional keyword args
	zero	= 'min'		if 'zero' not in kwargs else kwargs['zero']
	color	= 'purple'	if 'color' not in kwargs else kwargs['color']
	ax		= None		if 'ax' not in kwargs else kwargs['ax']

	### interpret input
	if isinstance(nbin, str) and nbin == 'auto':
		nbin = ars.optbins(A,50)
	elif not ars.isinteger(nbin):
		print("Error: Cannot calculate PMF - number of histogram bins must be either 'auto' or integer.\n")
		sys.exit()
	if isinstance(Alim_bin, str) and Alim_bin == 'auto':
		Alim_bin = [ min(A), max(A) ]
	elif not ars.isarray(Alim_bin) or len(Alim_bin) != 2:
		print("Error: Cannot calculate PMF - variable limits must be either 'auto' or 2-element array.\n")
		sys.exit()
	if isinstance(Alim_plot, str) and Alim_plot == 'auto':
		dAbin = (Alim_bin[1]-Alim_bin[0])/nbin
		Alim_plot = [ Alim_bin[0]-dAbin/2, Alim_bin[1]+dAbin/2 ]
	elif not ars.isarray(Alim_plot) or len(Alim_plot) != 2:
		if Alim_plot is not None:
			print("Error: Cannot calculate PMF - variable limits must be either 'auto' or 2-element array.\n")
			sys.exit()
	if zero != 'min' and zero != 'max':
		print("Error: Cannot calculate PMF - zero must be either 'min' or 'max'.\n")
		sys.exit()
	if ax is not None:
		if figLabel != 'auto' and figLabel is not None:
			print("Warining: Unused figure label input.")
			figLabel = None
	elif figLabel == 'auto':
		figLabel = "PMF"

	### calcualte PMF
	counts, bin_edges = np.histogram(A, nbin, range=Alim_bin)
	bins = np.ones(nbin)*np.nan
	PMF = np.ones(nbin)*np.nan
	for i in range(nbin):
		if counts[i] != 0:
			PMF[i] = -np.log(counts[i])
		bins[i] = (bin_edges[i]+bin_edges[i+1])/2
	if zero == 'min':
		PMF -= np.nanmin(PMF)
	elif zero == 'max':
		PMF -= np.nanmax(PMF)

	### initialize figure
	if ax is not None:
		plt.sca(ax)
	else:
		plt.figure(figLabel)

	### plot PMF
	plt.plot(bins, PMF, '-o', color=color)
	if Alim_plot is not None:
		plt.xlim(Alim_plot)
	if Alabel is not None:
		plt.xlabel(Alabel)
	plt.ylabel("PMF [kT]")
	if title is not None:
		plt.title(title)
	return bins, PMF


### plot umbrella sampling histograms and PMF
def plotUS(ops_eq, weights, ops_ts, ops_wham, PMF, PMF_err, opLabel="Order Parameter", titlePrefix=None, figLabelPrefix=None, nbin='auto', **kwargs):

	### additional keyword args
	insideLegend		= False		if 'insideLegend' not in kwargs else kwargs['insideLegend']
	op_precision		= 1			if 'op_precision' not in kwargs else kwargs['op_precision']
	weight_precision	= 1			if 'weight_precision' not in kwargs else kwargs['weight_precision']
	useOxUnits			= False		if 'useOxUnits' not in kwargs else kwargs['useOxUnits']
	setFigSize			= True		if 'setFigSize' not in kwargs else kwargs['setFigSize']

	### interpret input
	auto_bin = False
	if isinstance(nbin, str) and nbin == 'auto':
		auto_bin = True
	elif not ars.isinteger(nbin):
		print("Flag: Skipping US plot - number of histogram bins must be either 'auto' or integer.")
		return

	### adjust to oxDNA units (assuming nm and kcal)
	if useOxUnits:
		ops_eq *= 1/0.8518
		weights *= 1/8.2
		ops_ts = [i/0.8518 for i in ops_ts]
		ops_wham *= 1/0.8518

	### preparations
	nsim = len(ops_eq)
	titleHist = "US Histograms"
	titlePMF = "Free Energy"
	if titlePrefix is not None:
		titleHist = titlePrefix + " " + titleHist
		titlePMF = titlePrefix + " " + titlePMF
	figLabelHist = "US Hist"
	figLabelPMF = "PMF"
	if titlePrefix is not None:
		figLabelHist = f"{figLabelPrefix} {figLabelHist}"
		figLabelPMF = f"{figLabelPrefix} {figLabelPMF}"

	### plot nice histograms
	if setFigSize:
		plt.figure(figLabelHist, figsize=(10,6))
	else:
		plt.figure(figLabelHist)
	ax = plt.subplot(111)
	box = ax.get_position()
	ax.set_position([ box.x0, box.y0, box.width*0.8, box.height ])
	for i in range(nsim):
		if auto_bin == True:
			nbin = ars.optbins(ops_ts[i], 50)
		label = f"$OP={ops_eq[i]:0.{op_precision}f}$, $w={weights[i]:0.{weight_precision}f}$"
		ars.plotHist(ops_ts[i], figLabel=figLabelHist, nbin=nbin, Alim_plot=None, alpha=0.4, plotGauss=True, label=label)
	plt.xlabel(opLabel)
	plt.ylabel('Density')
	if insideLegend:
		plt.legend(loc='best')
	else:
		plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
	plt.title(titleHist)

	### plot PMF
	if setFigSize:
		plt.figure(figLabelPMF, figsize=(10,6))
	else:
		plt.figure(figLabelPMF)
	ax = plt.subplot(111)
	box = ax.get_position()
	ax.set_position([ box.x0, box.y0, box.width*0.8, box.height ])
	plt.errorbar(ops_wham, PMF, PMF_err, color='purple')
	plt.xlabel(opLabel)
	plt.ylabel('PMF [kcal/mol]')
	plt.title(titlePMF)


### calculate and plot mean squared displacement
def plotMSD(points, dbox3, dt_per_use, xLabel=None, yLabel=None, title=None, figLabel="MSD", nbin=10):

	### define center of mass function
	if len(points.shape) == 2:
		def getCOM(point):
			return point
	elif len(points.shape) == 3:
		def getCOM(points):
			return np.mean(points, axis=0)

	### calculate MSD
	nstep = int(points.shape[0])
	nstep_bin = int(np.floor(nstep/nbin))
	if nstep_bin < 2:
		print("Error: Cannot calculate MSD - too few steps per MSD bin.\n")
		sys.exit()
	time = np.linspace(0, (nstep_bin-1)*dt_per_use, nstep_bin)
	dis2 = np.zeros((nstep_bin, nbin))
	for q in range(nbin):
		com0 = getCOM(points[q*nstep_bin])
		com_prev = com0
		pbc = np.zeros(3)
		for i in range(1, nstep_bin):
			com = getCOM(points[q*nstep_bin+i])
			pbc += np.round((com - com_prev) / dbox3).astype(int)
			dis2[i,q] += np.sum((com - pbc*dbox3 - com0)**2)
			com_prev = com
	MSD = np.mean(dis2, axis=1)
	slope = linregress(time, MSD)[0]
	D = slope/6

	### plot MSD
	plt.figure(figLabel)
	plt.plot(time, MSD)
	plt.plot(time, slope*time)
	for q in range(min(nbin,10)):
		plt.plot(time, dis2[:,q], 'black', alpha=0.2)
	plt.xlabel(xLabel)
	plt.ylabel(yLabel)
	plt.legend(["Simulation","Linear Fit"], loc='upper left')
	plt.title(title)
	return D


### read file from lammps bond_write command, plot the energy and force
def plotBondWrite(bondWriteFile):

	### initialize
	lengths = []
	energies = []
	forces = []
		
	### open file, skip header
	ars.checkFileExist(bondWriteFile,"bond write")
	with open(bondWriteFile, 'r') as f:
		for line in range(5):
			next(f)
		
		### read the table
		for line in f:
			values = line.strip().split()
			if len(values) >= 2:
				lengths.append(float(values[1]))
				energies.append(float(values[2]))
				forces.append(float(values[3]))

	### plot setup
	fig, ax1 = plt.subplots(figsize=(8,6))
	ax2 = ax1.twinx()
	plt.title('Tabulated Bond')

	### energy plot
	color1 = 'tab:blue'
	ax1.set_xlabel('Length')
	ax1.set_ylabel('Energy', color=color1)
	ax1.plot(lengths, energies, color=color1)
	ax1.tick_params(axis='y', labelcolor=color1)

	### force plot
	color2 = 'tab:purple'
	ax2.set_ylabel('Force', color=color2)
	ax2.plot(lengths, forces, color=color2)
	ax2.tick_params(axis='y', labelcolor=color2)


################################################################################
### Calculations

### use wham to calculate PMF
def calcPMF(ops, weights, op_ts, tsFold, nbin, bin_padding, whamMetaFile, whamOutFile):

	### notes
	# the weights need to be in kcal/mol*{op unit}^2.
	# the output PMF are in kcal/mol.
	# uses autocorrelation to calculate the correlation time (used for error bars).

	### prepare timeseries file names
	ars.createEmptyFold(tsFold)

	### loop over simulations
	print("Writing timeseries...")
	nsim = len(ops)
	for i in range(nsim):
		nstep = op_ts[i].shape[0]

		### write sep timeseries file
		tsFile = tsFold + f"ts_sim{i:02.0f}"
		with open(tsFile, 'w') as f:
			for j in range(nstep):
				f.write(f"{j}\t{op_ts[i][j]}\n")

	### write wham metadata file
	print("Writing wham metadata...")
	with open(whamMetaFile, 'w') as f:
		for i in range(nsim):
			tsFile = tsFold + f"ts_sim{i:02.0f}"
			tau_int = 1 + 2*ars.calcCorrTime(op_ts[i])
			f.write(tsFile + f"\t{ops[i]}\t{weights[i]}\t{tau_int}\n")

	### wham parameters (not expected to change)
	T_wham = 300
	tol_wham = 0.001
	ntrial_MC_wham = 10

	### set wham bin limits
	bin_min = min(ops)-bin_padding
	bin_max = max(ops)+bin_padding

	### run wham (from Grossfield)
	print("Running wham...\n")
	subprocess.call(["wham", str(bin_min), str(bin_max), str(nbin), str(tol_wham), str(T_wham), "0", whamMetaFile, whamOutFile, str(ntrial_MC_wham), "37"])
	print("")

	### read wham output, return
	bins, PMF, PMF_err = ars.readWham(whamOutFile, nbin)
	return bins, PMF, PMF_err


### shift trajectory, placing the given point at the center, optionally unwrapping molecules at boundary
def centerPointsMolecule(points, molecules, dbox3s, center=1, unwrap=True, report_every=100):

	### notes
	# as the notation suggests, the indices contained in molecules must start at 1.
	# particles located perfectly at the origin are assumed to be dummy, and are thus
	  # not included in the center of mass calculation.

	### get counts
	nstep = points.shape[0]
	npoint = points.shape[1]
	nmolecule = int(max(molecules))

	### interpret input
	if ars.isnumber(dbox3s):
		dbox3s = np.ones(nstep)*dbox3s
	if ars.isarray(dbox3s) and ars.isnumber(dbox3s[0]) and len(dbox3s) == 3:
		dbox3s = np.ones((nstep,3))*dbox3s
	elif not ars.isarray(dbox3s) or len(dbox3s) != nstep or not (ars.isnumber(dbox3s[0]) or (ars.isarray(dbox3s[0]) and len(dbox3s[0]) != 3)):
		print("Error: Cannot center points - dbox3s must be number, nstep-element array, or nstep x 3 array.\n")
		sys.exit()

	### sort points by molecule
	points_moleculed = ars.sortPointsByMolecule(points, molecules)

	### initialize
	molecule_coms = np.zeros((nmolecule,3))
	points_centered = np.zeros((nstep,npoint,3))

	### loop over steps
	print("Centering trajectory...")
	for i in range(nstep):

		### calculate molecule coms
		for j in range(nmolecule):
			if ars.checkAllDummy(points_moleculed[j][i]):
				molecule_coms[j] = np.zeros(3)
			elif ars.checkAnyDummy(points_moleculed[j][i]):
				print("Error: Molecule contains mixed dummy and activated beads.\n")
				sys.exit()
			molecule_coms[j] = ars.calcCOM(points_moleculed[j][i], dbox3s[i])

		### set centering point
		if center == 'none':
			com = np.zeros(3)
		elif center == 'com' or center == 'com_points' or center == 'com_beads' or center == 'com_bases':
			com = ars.calcCOMnoDummy(points[i], dbox3s[i])
		elif center == 'com_molecules' or center == 'com_clusters':
			com = ars.calcCOMnoDummy(molecule_coms, dbox3s[i])
		elif ars.isinteger(center) and center <= nmolecule:
			com = molecule_coms[center-1]
		else:
			print("Error: Cannot center points - center must be either 'none', 'com', 'com_molecules', or integer <= nmolecule.\n")
			sys.exit()

		### center the points
		for j in range(npoint):
			if not all(molecule_coms[molecules[j]-1]==0):
				points_centered[i,j] = ars.applyPBC(points[i,j]-com, dbox3s[i])

		### unwrap molecules at boundary
		if unwrap:
			molecule_coms_centered = np.zeros((nmolecule,3))
			for j in range(nmolecule):
				if not all(molecule_coms[j]==0):
					molecule_coms_centered[j] = ars.applyPBC(molecule_coms[j]-com, dbox3s[i])
			for j in range(npoint):
				ref = molecule_coms_centered[molecules[j]-1]
				points_centered[i,j] = ref + ars.applyPBC(points_centered[i,j]-ref, dbox3s[i])

		### progress update
		if report_every and (i+1)%report_every == 0:
			print(f"centered {i+1} steps...")

	### result
	return points_centered


### shift trajectory, placing com of all beads at the center
def centerPointsBead(points, dbox3s):

	### notes
	# particles located perfectly at the origin are assumed to be dummy, and are thus
	  # not included in the center of mass calculation.

	### get counts
	nstep = points.shape[0]
	npoint = points.shape[1]

	### interpret input
	if ars.isnumber(dbox3s):
		dbox3s = np.ones(nstep)*dbox3s
	if ars.isarray(dbox3s) and ars.isnumber(dbox3s[0]) and len(dbox3s) == 3:
		dbox3s = np.ones((nstep,3))*dbox3s
	elif not ars.isarray(dbox3s) or len(dbox3s) != nstep or not (ars.isnumber(dbox3s[0]) or (ars.isarray(dbox3s[0]) and len(dbox3s[0]) != 3)):
		print("Error: Cannot center points - dbox3s must be number, nstep-element array, or nstep x 3 array.\n")
		sys.exit()

	### center points
	print("Centering trajectory...")
	points_centered = np.zeros((nstep,npoint,3))
	for i in range(nstep):
		com = ars.calcCOMnoDummy(points[i], dbox3s[i])
		for j in range(npoint):
			if not all(points[i,j]==0):
				points_centered[i,j] = ars.applyPBC(points[i,j]-com, dbox3s[i])
	return points_centered


### calculate center of mass, excluding dummy particle
def calcCOMnoDummy(r, dbox3):
	if ars.checkAllDummy(r):
		return np.zeros(3)
	mask_noDummy = ~np.all(r==0, axis=1)
	return ars.calcCOM(r[mask_noDummy], dbox3)


### calculate center of mass, using method from Bai and Breen 2008
def calcCOM(r, dbox3):
	xi_bar = np.mean( np.cos(2*np.pi*(r/dbox3+1/2)), axis=0 )
	zeta_bar = np.mean( np.sin(2*np.pi*(r/dbox3+1/2)), axis=0 )
	theta_bar = np.arctan2(-zeta_bar, -xi_bar) + np.pi
	r_ref = dbox3*(theta_bar/(2*np.pi)-1/2)
	com = r_ref + np.mean( ars.applyPBC(r-r_ref, dbox3), axis=0 )
	return com


### place each point in the same image as its preceeding neighbor
def unwrapChain(r, dbox3):
	npoint = r.shape[0]
	r_unwrapped = np.zeros((npoint,3))
	r_unwrapped[0] = r[0]
	for j in range(1,npoint):
		ref = r_unwrapped[j-1]
		r_unwrapped[j] = ref + ars.applyPBC( r[j]-ref, dbox3 )
	return r_unwrapped


### align the principal components of the given points with coordinate axes
def alignPCs(r, indices='all', axis_rank=[0,1,2], getPCs=False):

	### notes
	# indices gives the points from which to center, calculate PCs, and rotate.
	# axis rank gives the PC that each axis gets (in default, x-axis gets largest PC).
	# as is standard, the positions array contains points in the first dimension and
	  # coordinates in the second; in covariance terminology, this translates to 
	  # observations in rows and variables in columns.

	### interpret input
	if isinstance(indices, str) and indices == 'all':
		indices = np.arange(len(r))
	elif ars.isinteger(indices):
		indices = np.arange(indices)
	elif not ars.isarray(indices) or not ars.isinteger(indices[0]):
		print("Error: indices must be 'auto', integer, or int array.\n")
		sys.exit()

	### center about the given points
	com = np.mean(r[indices], axis=0)
	r_centered = r - com

	### get principal components
	cov = np.cov(r_centered[indices], rowvar=False)				
	eigenvalues, eigenvectors = np.linalg.eigh(cov)
	PCs_decreasing = eigenvectors[:,np.argsort(eigenvalues)[::-1]]
	PCs_axisRanked = PCs_decreasing[:,axis_rank]

	### enforce right-handedness
	if np.linalg.det(PCs_axisRanked) < 0:
		PCs_axisRanked[:, -1] *= -1

	### rotate positions, add back center
	r_rot = r_centered @ PCs_axisRanked
	r_aligned = r_rot + com

	### results
	output = [ r_aligned ]
	if getPCs: output.append(PCs_axisRanked)
	if len(output) == 1: output = output[0]
	return output


### histogram that blurs the lines between
def fractionalHist(data, bin_edges, alpha=1, limit_left='taper', limit_right='taper'):

	### notes
	# this script assumes constant bin width.
	# this method gives the same results as a normal histogram as alpha approaches infinity.
	# taper - treat like the other fractional bins, tapering bin_width/2 outside the limit.
	# cut - treat outside half like a canonical bin, with full attribution up to limit.
	# extend - full attribution from bin center past limit to infinity.

	### check input
	if alpha < 1:
		print("Error: alpha must be >= 1 for fractional histogram.\n")
		sys.exit()

	### metadata
	centers = (bin_edges[:-1] + bin_edges[1:]) / 2
	counts = np.zeros_like(centers)
	bin_width = bin_edges[1] - bin_edges[0]
	scale = 0.5**(alpha-1)

	### loop over data
	for x in data:

		# left limit
		if x <= centers[0]:
			if limit_left == 'taper':
				frac = (centers[0]-x) / bin_width
				if frac <= 0.5:
					counts[0] += 1 - (frac**alpha)/scale
				elif frac < 1:
					counts[0] += ((1-frac)**alpha)/scale
			elif limit_left == 'cut':
				if x > centers[0] - bin_width/2:
					counts[0] += 1.0
			elif limit_left == 'extend':
				counts[0] += 1.0
			else:
				print("Error: Unknown left end type, must be taper, cut, or extend.\n")
				sys.exit()
			continue

		# right limit
		if x >= centers[-1]:
			if limit_right == 'taper':
				frac = (x-centers[-1]) / bin_width
				if frac <= 0.5:
					counts[-1] += 1 - (frac**alpha)/scale
				elif frac < 1:
					counts[-1] += ((1-frac)**alpha)/scale
			elif limit_right == 'cut':
				if x < centers[-1] + bin_width/1:
					counts[-1] += 1.0
			elif limit_right == 'extend':
				counts[-1] += 1.0
			else:
				print("Error: Unknown right end type, must be taper, cut, or extend.\n")
				sys.exit()
			continue

		### regular bins
		for i in range(len(centers)-1):
			c_left = centers[i]
			c_right = centers[i+1]
			if c_left <= x <= c_right:
				frac_left = (x - c_left) / bin_width
				if frac_left <= 0.5:
					counts[i] += 1 - (frac_left**alpha)/scale
					counts[i+1] += (frac_left**alpha)/scale
				else:
					counts[i] += ((1-frac_left)**alpha)/scale
					counts[i+1] += 1 - ((1-frac_left)**alpha)/scale
				break

	### result
	return counts


### calculate optimum number of histogram bins
def optbins(A, maxM):
	N = len(A)
	logp = np.zeros(maxM)
	for M in range(1, maxM+1):
		n = np.histogram(A,bins=M)[0]
		part1 = N*np.log(M) + gammaln(M/2) - gammaln(N+M/2)
		part2 = -M*gammaln(1/2) + np.sum(gammaln(n+1/2))
		logp[M-1] = part1 + part2
	optM = np.argmax(logp) + 1
	return optM


### calculate sem using autocorrelation method
def calcSEMautocorr(A, hush=False):
	if np.isscalar(A):
		return 0
	tau = ars.calcCorrTime(A, hush)
	tau_int = 1 + 2*tau
	Neff = len(A)/tau_int
	return np.std(A) / np.sqrt(Neff)


### calculate correlation time
def calcCorrTime(A, hush=False):
	acf = ars.calcAutocorr(A)

	### Sokal window method
	M = 5
	tau = None
	tau_cumsum = np.cumsum(acf[1:])	
	for i, tau_est in enumerate(tau_cumsum):
		if i + 1 >= M * tau_est:
			tau = tau_est
			break

	### check for unconverged correlation
	if tau == None:
		print("Error: Correlation time estimate did not converge.")
		sys.exit()

	### check for uncorrelated data
	if tau < 1:
		tau = 0

	### check for poor estimate
	tau_int = 1+2*tau
	if not hush and len(A) < 50*tau_int:
		print("Warning: Correlation time estimate could be inaccurate.")
	
	### result
	return tau


### calculate autocorrelation function (slow but intuitive)
def calcAutocorrSlow(A):
	N = len(A)
	Ac = A - np.mean(A)
	acv = np.zeros(N)
	for l in range(N):
		acv[l] = np.dot(Ac[:N-l], Ac[l:N]) / (N-l)
	acf = acv / acv[0]
	return acf


### calculate autocorrelation function (fast but unintuitive)
def calcAutocorr(A):
	N = len(A)
	Ac = A - np.mean(A)
	fft_size = 2 ** int(np.ceil(np.log2(2 * N)))
	fft_data = np.fft.fft(Ac, n=fft_size)
	psd = fft_data * np.conj(fft_data)
	acv = np.fft.ifft(psd).real[:N]
	acv /= np.arange(N, 0, -1)
	acf = acv / acv[0]
	return acf


### calculate sem assuming independent measurements
def calcSEM(A):
	return np.std(A) / np.sqrt(len(A))


### check for overlap between a bead and a list of beads
def checkOverlap(r0, r_other, sigma, dbox3):
	return np.any(np.linalg.norm(ars.applyPBC(r_other-r0, dbox3),axis=1) < sigma)


### calculate moving average of an array
def movingAvg(A, stride=1):
	if stride < 1:
		print("Flag: Skipping moving average, stride must be positive integer.")
	avg = np.convolve(A, np.ones(stride)/stride, mode='valid')
	pad_left = (stride - 1) // 2
	pad_right = stride // 2
	return np.pad(avg, (pad_left,pad_right), mode='edge')



################################################################################
### Expressions

### divide vector by magnitude
def unitVector(vector):
	return vector / np.linalg.norm(vector)


### correlation between two vectors
def calcUnitDot(v0, v1):
	v0_u = ars.unitVector(v0)
	v1_u = ars.unitVector(v1)
	dot = np.clip(np.dot(v0_u, v1_u), -1, 1)
	return dot


### angle between two vectors
def calcAngle(v0, v1):
	return np.degrees(np.arccos(ars.calcUnitDot(v0, v1)))


### apply periodic boundary condition
def applyPBC(r, dbox):

	### notes
	# handles single values, single points, or arrays of points
	# dbox can be single value or ndim-element array
	# arrays must be numpy arrays

	return r - dbox*np.round(r/dbox)


### unit vector with random orientaiton and given number of dimensions
def randUnitVec(rng=np.random, ndim=3):
	return ars.unitVector(rng.normal(size=ndim))


### vector randomly placed within box
def randPos(dbox3, rng=np.random):
	if ars.isnumber(dbox3):
		dbox3 = np.ones(3)*dbox3
	x = np.zeros(3)
	for i in range(3):
		x[i] = rng.uniform(-dbox3[i]/2, dbox3[i]/2)
	return x


### calculate weighted mean
def calcMean(A, weights=None):
	if weights is None:
		return np.mean(A)
	else:
		return np.sum(weights*A)/np.sum(weights)


### calculate weighted standard deviation
def calcNormStats(A, weights=None):
	if weights is None:
		return np.mean(A), np.std(A)
	else:
		sumw = np.sum(weights)
		mu = np.sum(weights*A)/sumw
		sigma = np.sqrt(np.sum(weights*(A-mu)**2)/sumw)
		return mu, sigma


################################################################################
### Random

### get common nice colors
def getColor(color):
	if color == 'teal':
		return np.array([0,145,147])
	elif color == 'orchid':
		return np.array([122,129,255])
	elif color == 'silver':
		return np.array([214,241,241])
	elif color == 'purple':
		return np.array([68,1,84])
	elif color == 'grey':
		return np.array([153,153,153])
	else:
		print("Error: Unknown color.\n")
		sys.exit()


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


### determine if position array contains any dummy beads
def checkAllDummy(r):
	return np.all(r==0)


### trim op and PMF array to match finite values of PMF array
def trimUS(op, PMF, PMF_err):
	finite_idxs = np.isfinite(PMF)
	start_idx = np.argmax(finite_idxs)
	end_idx = len(PMF) - np.argmax(finite_idxs[::-1])
	op_trimmed = op[start_idx:end_idx]
	PMF_trimmed = PMF[start_idx:end_idx]
	PMF_err_trimmed = PMF_err[start_idx:end_idx]
	return op_trimmed, PMF_trimmed, PMF_err_trimmed


### determine if file exists
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


### creates new folder, only if it doesn't already exist
def createSafeFold(newFold):
	os.makedirs(newFold, exist_ok=True)


### creates new empty folder
def createEmptyFold(newFold):
	if os.path.exists(newFold):
		shutil.rmtree(newFold)
	os.makedirs(newFold)


### render sequential array compactly
def compressSeqArr(A):
	if all(A[i]+1 == A[i+1] for i in range(len(A)-1)):
		return f"[ {A[0]} ... {A[-1]} ]"
	else:
		return str(A)


### test if variable is numeric (both float and integer count)
def isnumber(x):
	try:
		value = float(x)
		return True
	except:
		return False


### test if variable is an integer (both python int and numpy int count)
def isinteger(x):
	if isinstance(x, int):
		return True
	elif isinstance(x, np.int64):
		return True
	else:
		return False


### check if variable is an array (both list and numpy array count)
def isarray(x):
	if isinstance(x, list):
		return True
	elif isinstance(x, np.ndarray):
		return True
	else:
		return False




