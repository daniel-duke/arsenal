import numpy as np
import sys
import os

## Description
# this script reads a lammps-style atom dump, centers the atoms (either
  # by atom com, molecule com, or given molecule ID), optionally unwraps
  # the atoms from the boundary by molecule, and optionally sets the color
  # in the new trajectory according to the molecule ID.
# "col2s" refers to whatever is in the second column of the trajectory file;
  # for a regular atom dump, this is the atom type, but it could be anything
  # for a custom dump.
# this script runs from the command line with four argumnets: (1) the name of
  # the trajectory file to load and center, (2) what to center, (3) whether
  # to unwrap the trajecotry, and (4) whether to color the trajectory.
# the script outputs a file with the same name as the input file, modified to
  # with "_centered" at the end (but before ".dat"); for example, the input file
  # "trajectory.dat" would result in output file "trajectory_centered.dat"
# this script assumes the input particle positions are written in the "scaled"
  # format, meaning the actual positions equal the value times the box diameter.


################################################################################
### Parameters

def main():

    ### get command line input
    datFile = sys.argv[1]
    center = sys.argv[2]
    unwrap = int(sys.argv[3])
    set_color = int(sys.argv[4])

    ### optional fifth argument
    display = False
    if len(sys.argv) > 5:
        display = int(sys.argv[5])

    ### read file (no skipped steps, no coarsening)
    points, col2s, dbox3s, dump_every = readAtomDump(datFile, display)

    ### center the points
    points_centered = centerPointsMolecule(points, col2s, dbox3s, center, unwrap, display)

    ### write output file
    outDatFile = datFile[:(len(datFile)-4)] + "_centered" + datFile[(len(datFile)-4):]
    writeAtomDump(outDatFile, points_centered, col2s, dbox3s, dump_every, set_color)


################################################################################
### File Managers

### read lammps-style trajectory
def readAtomDump(datFile, display):
    if display: print("Loading LAMMPS-style trajectory...")
    testFileExist(datFile, "trajectory")
    with open(datFile, 'r') as f:
        content = f.readlines()
    if display: print("Parsing trajectory...")

    nbd_total = int(content[3].split()[0])
    dump_every = int(content[nbd_total+10].split()[0]) - int(content[1].split()[0])
    nstep = int(len(content)/(nbd_total+9))

    points = np.zeros((nstep,nbd_total,3))
    col2s = np.zeros(nbd_total,dtype=int)
    dbox3s = np.zeros((nstep,3))
    for i in range(nstep):
        for k in range(3):
            line = content[(nbd_total+9)*i+5+k].split()
            dbox3s[i,k] = float(line[1]) - float(line[0])
        for j in range(nbd_total):
            line = content[(nbd_total+9)*i+9+j].split()
            if i == 0:
                col2s[j] = int(line[1])
            points[i,j] = (np.array(line[2:5],dtype=float)-1/2)*dbox3s[i]
        points[i] = applyPBC(points[i], dbox3s[i])
        if display and i%1000 == 0 and i != 0:
            print(f"Processed {i} steps...")
    return points, col2s, dbox3s, dump_every


### write lammps-style atom dump
def writeAtomDump(outDatFile, points, col2s, dbox3s, dump_every, set_color):
    nstep = points.shape[0]
    npoint = points.shape[1]
    len_npoint = len(str(npoint))
    len_ncol2 = len(str(max(col2s)))
    with open(outDatFile,'w') as f:
        for i in range(nstep):
            len_dbox = len(str(int(max(dbox3s[i]))))
            f.write(f"ITEM: TIMESTEP\n{i*dump_every}\n")
            f.write(f"ITEM: NUMBER OF ATOMS\n{npoint}\n")
            f.write(f"ITEM: BOX BOUNDS pp pp pp\n")
            f.write(f"-{dbox3s[i,0]/2:0{len_dbox+3}.2f} {dbox3s[i,0]/2:0{len_dbox+3}.2f} xlo xhi\n")
            f.write(f"-{dbox3s[i,1]/2:0{len_dbox+3}.2f} {dbox3s[i,1]/2:0{len_dbox+3}.2f} ylo yhi\n")
            f.write(f"-{dbox3s[i,2]/2:0{len_dbox+3}.2f} {dbox3s[i,2]/2:0{len_dbox+3}.2f} zlo zhi\n")
            if set_color:
                f.write("ITEM: ATOMS id type xs ys zs\n")
            else:
                f.write("ITEM: ATOMS id mol xs ys zs\n")
            for j in range(npoint):
                f.write(f"{j+1:<{len_npoint}} " + \
                        f"{col2s[j]:<{len_ncol2}}  " + \
                        f"{points[i,j,0]/dbox3s[i,0]+1/2:10.8f} " + \
                        f"{points[i,j,1]/dbox3s[i,1]+1/2:10.8f} " + \
                        f"{points[i,j,2]/dbox3s[i,2]+1/2:10.8f}\n")


################################################################################
### Calculations

### shift trajectory, placing the given point at the center, optionally unwrapping molecules at boundary
def centerPointsMolecule(points, molecules, dbox3s, center, unwrap, display):
    nstep = points.shape[0]
    npoint = points.shape[1]
    nmolecule = int(max(molecules))

    ### sort points by molecule
    points_moleculed = sortPointsByMolecule(points, molecules)

    ### initialize
    molecule_coms = np.zeros((nmolecule,3))
    points_centered = np.zeros((nstep,npoint,3))

    ### loop over steps
    if display: print("Centering trajectory...")
    for i in range(nstep):

        ### calculate molecule coms
        for j in range(nmolecule):
            if checkAllDummy(points_moleculed[j][i]):
                molecule_coms[j] = np.zeros(3)
            elif checkAnyDummy(points_moleculed[j][i]):
                print("Error: Molecule contains mixed dummy and activated beads.\n")
                sys.exit()
            molecule_coms[j] = calcCOM(points_moleculed[j][i], dbox3s[i])

        ### set centering point
        if center == 'none':
            com = np.zeros(3)
        elif center == 'com_points' or center == 'com_beads' or center == 'com_bases':
            com = calcCOMnoDummy(points[i], dbox3s[i])
        elif center == 'com_molecules' or center == 'com_clusters':
            com = calcCOMnoDummy(molecule_coms, dbox3s[i])
        elif isnumber(center) and int(center) <= nmolecule:
            com = molecule_coms[int(center)-1,:]
        else:
            print("Error: Cannot center points - center must be either 'none', 'com_points', 'com_molecules', or integer <= nmolecule.\n")
            sys.exit()

        ### center the points
        for j in range(npoint):
            if not all(molecule_coms[molecules[j]-1]==0):
                points_centered[i,j] = applyPBC(points[i,j]-com, dbox3s[i])

        ### unwrap molecules at boundary
        if unwrap:
            molecule_coms_centered = np.zeros((nmolecule,3))
            for j in range(nmolecule):
                if not all(molecule_coms[j]==0):
                    molecule_coms_centered[j] = applyPBC(molecule_coms[j]-com, dbox3s[i])
            for j in range(npoint):
                ref = molecule_coms_centered[molecules[j]-1]
                points_centered[i,j] = ref + applyPBC(points_centered[i,j]-ref, dbox3s[i])

        ### progress update
        if display and (i+1)%1000 == 0:
            print(f"Centered {i+1} steps...")

    ### result
    return points_centered


### calculate center of mass, excluding dummy particle
def calcCOMnoDummy(r, dbox3):
    if checkAllDummy(r):
        return np.zeros(3)
    mask_noDummy = ~np.all(r==0, axis=1)
    return calcCOM(r[mask_noDummy], dbox3)


### calculate center of mass, using method from Bai and Breen 2008
def calcCOM(r, dbox3):
    xi_bar = np.mean( np.cos(2*np.pi*(r/dbox3+1/2)), axis=0 )
    zeta_bar = np.mean( np.sin(2*np.pi*(r/dbox3+1/2)), axis=0 )
    theta_bar = np.arctan2(-zeta_bar, -xi_bar) + np.pi
    r_ref = dbox3*(theta_bar/(2*np.pi)-1/2)
    com = r_ref + np.mean( applyPBC(r-r_ref, dbox3), axis=0 )
    return com


################################################################################
### Calculations

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
            points_moleculed[molecules[j]-1][i,n_molecule_count[molecules[j]-1],:] = points[i,j,:]
            n_molecule_count[molecules[j]-1] += 1
    return points_moleculed


### determine if position array contains any dummy beads
def checkAnyDummy(r):
    return np.any(np.all(r==0, axis=1))


### determine if position array contains any dummy beads
def checkAllDummy(r):
    return np.all(r==0)


### test if files exist
def testFileExist(file, name):
    if not os.path.isfile(file):
        print("Error: Could not find " + name + " file:")
        print(file + "\n")
        sys.exit()


### test if variable is numeric (both float and integer count)
def isnumber(x):
    try:
        value = float(x)
        return True
    except:
        return False


### apply periodic boundary condition
def applyPBC(r, dbox):
    return r - dbox*np.round(r/dbox)


### run the script
if __name__ == "__main__":
    main()

