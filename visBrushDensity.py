import arsenal as ars
import numpy as np
from scipy.ndimage import gaussian_filter
import argparse
import sys

## Description
# this script reads an oxDNA trajectory of a nanostructure with floppy ssDNA
  # brushes and computes a 3D density map of where the brush nucleotides go
  # over the course of the trajectory.
# the cluster file contains the base indices for the brushes (one cluster per
  # brush); all of the remaining base indices make up the core.
# every frame is aligned to the mean structure, as controlled by three options:
  # --fit gives the nucleotides that the rotation is fit to (core, all, or
  # brushes), --axis restricts the rotation (free, none, or about a single
  # axis), and --pivot gives the point the rotation is about (by default the
  # center of mass of the fit nucleotides, which also aligns translation; or a
  # fixed point, e.g. for a structure that is pinned in place). a group of
  # nucleotides only has to exist if the alignment actually uses it, so there is
  # no need for a core when there is no rotation (--axis 0) about a fixed point
  # (--pivot 0).
# the density grid(s) are written as Gaussian cube files, which can be rendered
  # in ovito using the Create Isosurface modifier.
# the reference structure is written separately as a LAMMPS-style geometry file,
  # with the option to either represent nucleotides explicitly (default) or as
  # simple beads at their COMs.


################################################################################
### Parameters

def main():

	### get arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('topFile',			type=str, default=None,		help='name of topology file')
	parser.add_argument('datFile',			type=str, default=None,		help='name of configuration/trajectory file')
	parser.add_argument('clusterFile',		type=str, default=None,		help="name of cluster file")
	parser.add_argument('--nstep-skip', 	type=int, default=0,		help="number of recorded initial steps to skip")
	parser.add_argument('--coarse-time',	type=int, default=1,		help="coarse factor for time steps")
	parser.add_argument('--voxel',			type=float, default=0.5,	help="grid spacing, in nm")
	parser.add_argument('--padding',		type=float, default=2.0,	help="padding around the brush bounding box, in nm")
	parser.add_argument('--sigma',			type=float, default=0.5,	help="Gaussian smoothing width, in nm (0 to disable)")
	parser.add_argument('--fit',			type=str, default='core', choices=['core','all','brushes'], help="nucleotides that the rotation is fit to when aligning the frames")
	parser.add_argument('--axis',			type=str, nargs='+', default=None, help="rotation allowed when aligning the frames: leave out for free rotation, 0 for no rotation, 1/2/3 to only rotate about x/y/z, or three numbers to only rotate about that axis")
	parser.add_argument('--pivot',			type=str, nargs='+', default=None, help="what the frames are rotated about when aligned: com-core, com-all, or com-brushes rotate about the center of mass of that group and align translation using it; a number (0 for the origin) or three numbers X Y Z (in nm) rotate about that point and apply no translation (e.g. a pinned particle); by default the center of mass of the fit nucleotides")
	parser.add_argument('--alignToBox',		action='store_true',		help="rotate the output (for visualization only) so the principal components of the fit nucleotides' mean structure line up with the box axes")
	parser.add_argument('--noSplit',		action='store_true',		help="plot nucleotides as single bead at COM")
	parser.add_argument('--writeOvito',		action='store_true',		help="write an ovito scene file")
	args = parser.parse_args()

	### interpret the alignment options
	axis = None if args.axis is None else getArgType(args.axis)
	pivot = f"com-{args.fit}" if args.pivot is None else getArgType(args.pivot)
	if ars.isnumber(pivot): pivot = np.full(3, float(pivot))


################################################################################
### Heart

	### output folder
	outFold = "vis_brushes/"
	outGeoFile = outFold + "geometry.in"
	outOvitoFile = outFold + "brushes.ovito"

	### read topology and cluster definitions
	strands, bonds, nba_total = readTop(args.topFile)
	clusters = ars.readCluster(args.clusterFile)

	### identify the groups of nucleotides, and the ones the alignment uses
	groups = getGroups(clusters, nba_total)
	bais_fit, bais_pivot = getAlignmentGroups(groups, args.fit, axis, pivot, args.alignToBox)
	bais_brush = groups['brushes']

	### read trajectory, fixing periodic boundary breaks (and centering, unless the alignment pivots about a fixed point)
	points, dbox3, axes = ars.readOxDNA(args.datFile, args.nstep_skip, args.coarse_time, getDbox3=True, getAxes=True)
	points = ars.centerPointsMolecule(points, strands, dbox3, center='none' if bais_pivot is None else 'com', unwrap=True)

	### align every frame to the mean structure
	points, pivot = setPivot(points, bais_pivot, pivot)
	points_aligned, r_ref, axes_ref = alignFrames(points, axes, bais_fit, axis, pivot)

	### for visualization, rotate everything so the principal components of the fit nucleotides' mean structure line up with the box
	if args.alignToBox:
		pivot_box = np.mean(r_ref[bais_fit], axis=0)
		r_ref, PCs = ars.alignPCs(r_ref, bais_fit, getPCs=True)
		points_aligned = (points_aligned - pivot_box) @ PCs + pivot_box
		axes_ref = axes_ref @ PCs

	### assign colors from cluster identity (1 for the core, 2, 3, ... for each brush)
	colors = np.ones(nba_total, dtype=int)
	for c,cluster in enumerate(clusters):
		colors[np.array(cluster)] = c+2

	### create output folder (only now, so that a bad input does not wipe the previous results)
	ars.createEmptyFold(outFold)

	### write the reference structure: the mean shape, with backbone bonds
	if args.noSplit:
		ars.writeGeo(outGeoFile, dbox3, r_ref, strands, colors, bonds)

	else:
		### write geometry
		r_ref, quats_ref, radii, bonds, [colors] = ars.splitNucleotides(r_ref, axes_ref, bonds, [colors])
		ars.writeGeo(outGeoFile, dbox3, r_ref, types=colors, bonds=bonds, radii=radii, quats=quats_ref)

	### write the density grids
	outCubeFiles = writeDensities(outFold, points_aligned, bais_brush, clusters, args.voxel, args.padding, args.sigma)

	### write ovito scene
	if args.writeOvito:
		writeOvito(outOvitoFile, outGeoFile, outCubeFiles)


################################################################################
### File Handlers

### read oxdna topology
def readTop(topFile):
	ars.checkFileExist(topFile, "topology")
	with open(topFile) as f:
		content = f.readlines()
	nba_total = int(content[0].split()[0])
	strands = np.zeros(nba_total,dtype=int)
	bonds = np.ones((nba_total,3),dtype=int)
	bond_count = 0
	for i in range(nba_total):
		line = content[i+1].split()
		strands[i] = int(line[0])
		if int(content[i+1].split()[3]) != -1:
			bonds[bond_count,1] = i+1
			bonds[bond_count,2] = int(line[3])+1
			bond_count += 1
	bonds = bonds[:bond_count,:]
	return strands, bonds, nba_total


### write Gaussian cube file (for ovito)
def writeCube(cubeFile, grid, origin, spacing):

	### notes
	# a single placeholder atom is always written at the grid origin
	  # because ovito refuses to recognize the format of a cube file 
	  # with no atoms.
	# ovito assumes cube files use bohr length units and thus
	  # automatically converts the units to Angstrom; although this
	  # conversion can be deselected for the density values, the same
	  # is not true for the origin; thus, all length values written
	  # to the cube file are scaled by the Angstrom-to-bohr conversion
	  # factor, so that they may all be rescaled back to nm by ovito.

	### constants
	ang2bohr = 1.8897

	### counts
	nx, ny, nz = grid.shape

	### write file
	with open(cubeFile, 'w') as f:

		### header
		f.write("Brush density data\n\n")
		ox, oy, oz = np.asarray(origin)*ang2bohr
		f.write(f"{1:5d} {ox:12.6f} {oy:12.6f} {oz:12.6f}\n")
		dx, dy, dz = np.asarray(spacing)*ang2bohr
		f.write(f"{nx:5d} {dx:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
		f.write(f"{ny:5d} {0.0:12.6f} {dy:12.6f} {0.0:12.6f}\n")
		f.write(f"{nz:5d} {0.0:12.6f} {0.0:12.6f} {dz:12.6f}\n")

		### placeholder atom
		f.write(f"{0:5d} {0.0:12.6f} {ox:12.6f} {oy:12.6f} {oz:12.6f}\n")

		### density data
		for ix in range(nx):
			for iy in range(ny):
				line_vals = grid[ix,iy,:]
				for k in range(0,nz,6):
					chunk = line_vals[k:k+6]
					f.write(" ".join(f"{v:13.5e}" for v in chunk) + "\n")


### write the density grid of all the brushes together, then of each brush individually, all on the same grid
def writeDensities(outFold, points_aligned, bais_brush, clusters, voxel, padding, sigma):
	nstep = points_aligned.shape[0]
	nbrush = len(clusters)
	ars.initStatusBar("Writing cubes")

	### pool brush positions across all frames and determine the grid
	points_brush = points_aligned[:,bais_brush].reshape(-1,3)
	edges = calcGridEdges(points_brush, voxel, padding)

	### all brushes together
	grid, origin = calcDensityGrid(points_brush, nstep, edges, voxel, sigma)
	writeCube(outFold + "brushes.cube", grid, origin, np.full(3,voxel))
	ars.updateStatusBar(0,nbrush+1,units="brushes")

	### each brush individually
	outCubeFiles = []
	for c,cluster in enumerate(clusters):
		points_brush = points_aligned[:,combineClusters([cluster])].reshape(-1,3)
		grid, origin = calcDensityGrid(points_brush, nstep, edges, voxel, sigma)
		outCubeFiles.append(outFold + f"brush{c+1}.cube")
		writeCube(outCubeFiles[c], grid, origin, np.full(3,voxel))
		ars.updateStatusBar(c+1,nbrush+1,units="brushes")

	return outCubeFiles


### build an ovito scene with all pipelies already configured
def writeOvito(sceneFile, geoFile, cubeFiles):

	### import ovito
	from ovito.io import import_file
	from ovito.modifiers import CreateIsosurfaceModifier
	from ovito.modifiers import EditTypesModifier
	from ovito import scene

	### brush color sequence
	isoColors = [ ars.getColor('sky', scaled=True),
				  ars.getColor('orchid', scaled=True),
				  ars.getColor('teal', scaled=True),
				  ars.getColor('purple', scaled=True) ]

	### reference structure
	pipeline = import_file(geoFile)
	pipeline.source.data.cell.vis.enabled = False
	pipeline.source.data.particles.vis.radius = 0.34
	pipeline.source.data.particles.bonds.vis.width = 0.17

	### reference structure colors
	mod = EditTypesModifier()
	pipeline.modifiers.append(mod)
	mod.edit_type(1).color = ars.getColor('mercury',scaled=True)
	for t in range(1,len(cubeFiles)+1):
		mod.edit_type(t+1).color = ars.getColor('steel',scaled=True)
	pipeline.add_to_scene()

	### loop over brushes
	for i,cubeFile in enumerate(cubeFiles):

		### density grid
		pipeline = import_file(cubeFile)
		pipeline.source.data.cell.vis.enabled = False
		pipeline.source.data.particles.vis.enabled = False

		### isosurface
		mod = CreateIsosurfaceModifier(operate_on='voxels:imported', property='Property', isolevel=0.1, smoothing_level=10)
		pipeline.modifiers.append(mod)
		mod.vis.surface_color = isoColors[i%len(isoColors)]
		mod.vis.surface_transparency = 0.5
		pipeline.add_to_scene()

	### save file
	scene.save(sceneFile)


################################################################################
### Calculations

### determine voxel grid edges that bound the given points, with padding
def calcGridEdges(points, voxel, padding):
	mins = points.min(axis=0) - padding
	maxs = points.max(axis=0) + padding
	nvoxel = np.maximum(1, np.ceil((maxs-mins)/voxel).astype(int))
	return [ np.linspace(mins[d], mins[d]+nvoxel[d]*voxel, nvoxel[d]+1) for d in range(3) ]


### build a smoothed, normalized 3D density histogram from a point cloud
def calcDensityGrid(points, nstep, edges, voxel, sigma):

	### histogram the points
	grid = np.histogramdd(points, bins=edges)[0]

	### smooth
	if sigma > 0:
		grid = gaussian_filter(grid, sigma=sigma/voxel)

	### normalize to a number density (per frame, per nm^3)
	grid /= nstep*voxel**3

	### result
	origin = np.array([ e[0] for e in edges ])
	return grid, origin


### align every frame to the mean structure, averaging the orientations along the way
def alignFrames(points, axes, bais_fit, axis, pivot):
	nstep, nba_total = points.shape[:2]

	### calculate the mean structure, then align every frame to it
	r_ref = ars.calcMeanStructure(points, bais_fit, axis=axis, pivot=pivot)
	ars.initStatusBar("Aligning")
	points_aligned = np.zeros((nstep,nba_total,3))
	a1_sum = np.zeros((nba_total,3))
	a3_sum = np.zeros((nba_total,3))
	for i in range(nstep):
		points_aligned[i], R = ars.kabschAlgorithm(points[i], r_ref, indices=bais_fit, getR=True, axis=axis, pivot=pivot)
		a1_sum += axes[i,:,0] @ R
		a3_sum += axes[i,:,2] @ R
		ars.updateStatusBar(i,nstep)

	### renormalize the averaged orientations
	a1_ref = a1_sum / np.linalg.norm(a1_sum, axis=1, keepdims=True)
	a3_ref = a3_sum / np.linalg.norm(a3_sum, axis=1, keepdims=True)
	a2_ref = np.cross(a3_ref, a1_ref)
	a2_ref /= np.linalg.norm(a2_ref, axis=1, keepdims=True)
	axes_ref = np.stack([a1_ref, a2_ref, a3_ref], axis=1)

	### results
	return points_aligned, r_ref, axes_ref


################################################################################
### Utility Functions

### interpret a command line argument given as a keyword, a number, or a list of numbers
def getArgType(arg):
	arg = [ int(a) if a.isdigit() else float(a) if ars.isnumber(a) else a for a in arg ]
	return arg[0] if len(arg) == 1 else arg


### identify the groups of nucleotides (as lists of base indices): the brushes are all the clusters, and the core is everything else
def getGroups(clusters, nba_total):
	if len(clusters) == 0:
		print("Error: No brush clusters found in the cluster file.")
		sys.exit(1)
	bais_all = np.arange(nba_total)
	bais_brush = combineClusters(clusters)
	return { 'core'		: np.setdiff1d(bais_all, bais_brush),
			 'all'		: bais_all,
			 'brushes'	: bais_brush }


### get the groups of nucleotides that the alignment uses, checking that they exist when they are needed:
  # the fit group is needed to fit the rotation (or to find the principal components for aligning to
  # the box), and the pivot group is needed when the pivot is the center of mass of a group
def getAlignmentGroups(groups, fit, axis, pivot, alignToBox):

	### fit group (the rotation is not fit if there is no rotation, so anything can stand in)
	bais_fit = groups['all']
	noRotation = ars.isnumber(axis) and float(axis) == 0
	if not noRotation or alignToBox:
		bais_fit = groups[fit]
		if len(bais_fit) == 0:
			print(f"Error: There are no {fit} nucleotides, which are needed to fit the alignment; choose another group with --fit, or avoid fitting with --axis 0 (and no --alignToBox).")
			sys.exit(1)

	### pivot group (none if the pivot is a fixed point)
	bais_pivot = None
	if isinstance(pivot, str):
		if not pivot.startswith('com-') or pivot[4:] not in groups:
			print("Error: Unrecognized pivot, keywords are " + ", ".join("com-"+g for g in groups) + ".")
			sys.exit(1)
		bais_pivot = groups[pivot[4:]]
		if len(bais_pivot) == 0:
			print(f"Error: There are no {pivot[4:]} nucleotides, which are needed to find the pivot; choose another group or a fixed point with --pivot.")
			sys.exit(1)

	### results
	return bais_fit, bais_pivot


### if the pivot is the center of mass of a group, put that at the origin in every frame and pivot about the origin
  # (this aligns translation along with rotation)
def setPivot(points, bais_pivot, pivot):
	if bais_pivot is None:
		return points, pivot
	return points - np.mean(points[:,bais_pivot], axis=1, keepdims=True), np.zeros(3)


### get a single list of all base indices included in a list of clusters
def combineClusters(clusters):
	indices = { i for cluster in clusters for i in cluster }
	return np.array(sorted(set(indices)))


### run the script
if __name__ == "__main__":
	main()
	print()

