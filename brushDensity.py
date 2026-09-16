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
  # brush), with the the core structure assumed to comprise of all the non-brush
  # base indices.
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
	parser.add_argument('--noSplit',		action='store_true',		help="plot nucleotides as single bead at COM")
	parser.add_argument('--align',			action='store_true',		help="rotate mean structure to align its principal components with the box")
	parser.add_argument('--writeOvito',		action='store_true',		help="write an ovito scene file")
	args = parser.parse_args()


################################################################################
### Heart

	### output folder
	outFold = "vis_brushes/"
	outGeoFile = outFold + "geometry.in"
	outOvitoFile = outFold + "brushes.ovito"

	### create output folder
	ars.createEmptyFold(outFold)

	### read topology and cluster definitions
	strands, bonds, nba_total = readTop(args.topFile)
	clusters = ars.readCluster(args.clusterFile)
	bais_brush = combineClusters(clusters)
	bais_all = np.arange(nba_total)
	bais_core = np.setdiff1d(bais_all, bais_brush)
	nbrush = len(clusters)

	### read and center trajectory
	points, dbox3, axes = ars.readOxDNA(args.datFile, args.nstep_skip, args.coarse_time, getDbox3=True, getAxes=True)
	points = ars.centerPointsMolecule(points, strands, dbox3, center='com', unwrap=True)
	nstep = points.shape[0]

	### calculate reference structure
	r_ref = ars.calcMeanStructure(points, bais_core)
	if args.align: r_ref = ars.alignPCs(r_ref, bais_core)

	### calculate mean structure
	ars.initStatusBar("Aligning")
	points_aligned = np.zeros((nstep,nba_total,3))
	a1_sum = np.zeros((nba_total,3))
	a3_sum = np.zeros((nba_total,3))
	for i in range(nstep):
		points_aligned[i], R = ars.kabschAlgorithm(points[i], r_ref, indices=bais_core, getR=True)
		a1_sum += axes[i,:,0] @ R
		a3_sum += axes[i,:,2] @ R
		ars.updateStatusBar(i,nstep)

	### renormalize the averaged orientations
	a1_ref = a1_sum / np.linalg.norm(a1_sum, axis=1, keepdims=True)
	a3_ref = a3_sum / np.linalg.norm(a3_sum, axis=1, keepdims=True)
	a2_ref = np.cross(a3_ref, a1_ref)
	a2_ref /= np.linalg.norm(a2_ref, axis=1, keepdims=True)
	axes_ref = np.stack([a1_ref, a2_ref, a3_ref], axis=1)

	### assign colors from cluster identity
	colors = np.ones(nba_total, dtype=int)
	for c,cluster in enumerate(clusters):
		colors[np.array(cluster)] = c+2

	### write the reference structure: the core's mean shape, with backbone bonds
	if args.noSplit:
		ars.writeGeo(outGeoFile, dbox3, r_ref, strands, colors, bonds)

	else:
		### write geometry
		r_ref, quats_ref, colors, radii, bonds = splitNucleotides(r_ref, axes_ref, colors, bonds)
		ars.writeGeo(outGeoFile, dbox3, r_ref, types=colors, bonds=bonds, radii=radii, quats=quats_ref)

	### cube updates
	ars.initStatusBar("Writing cubes")

	### pool brush positions across all frames and determine grid
	points_brush = points_aligned[:,bais_brush].reshape(-1,3)

	### determine grid
	edges = calcGridEdges(points_brush, args.voxel, args.padding)

	### calculate density grid
	grid, origin = calcDensityGrid(points_brush, nstep, edges, args.voxel, args.sigma)

	### write cube file
	outCubeFile = outFold + "brushes.cube"
	writeCube(outCubeFile, grid, origin, np.full(3,args.voxel))
	ars.updateStatusBar(0,nbrush+1,units="brushes")

	### repeat for each brush individually
	outCubeFiles = []
	for c,cluster in enumerate(clusters):
		bais_brush = np.array(sorted(set(cluster)))

		### pool brushes, calculate grid
		points_brush = points_aligned[:,bais_brush].reshape(-1,3)
		grid, origin = calcDensityGrid(points_brush, nstep, edges, args.voxel, args.sigma)

		### write cube file
		outCubeFiles.append(outFold + f"brush{c+1}.cube")
		writeCube(outCubeFiles[c], grid, origin, np.full(3,args.voxel))
		ars.updateStatusBar(c+1,nbrush+1)

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


################################################################################
### Utility Functions

### get a single list of all base indices included in a list of clusters
def combineClusters(clusters):
	indices = { i for cluster in clusters for i in cluster }
	return np.array(sorted(set(indices)))


### split particles into backbone and base sites, adding appropriate bonds
def splitNucleotides(r, axes, colors, bonds):
	npoint = r.shape[0]

	### initialize
	r_split = np.zeros((npoint*2,3))
	quats = np.zeros((npoint*2,4))
	quats[:,0] = 1

	### split position data
	r_split[:npoint] = r + ars.ox2nm*(-0.34*axes[:,0] + 0.3408*axes[:,1])
	r_split[npoint:] = r + ars.ox2nm*(0.4*axes[:,0])

	### calculate orientation data
	for j in range(npoint):
		quats[npoint+j] = ars.axesToQuat(axes[j])

	### add base to backbone bonds
	for j in range(npoint):
		bonds = np.append(bonds,[[2,j+1,npoint+j+1]], axis=0)

	### identification data
	colors = np.concatenate((colors,colors))

	### bead sizes
	radii = np.full((npoint*2,3), 0.34)
	radii[npoint:,2] = 0.17

	### results
	return r_split, quats, colors, radii, bonds


### run the script
if __name__ == "__main__":
	main()
	print()

