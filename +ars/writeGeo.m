%%% write lammps geometry file
function writeGeo(geoFile,dbox,atoms,bonds,angles,options)

    %%% set arguments
    arguments
        geoFile string
        dbox double
        atoms double
        bonds double
        angles double
        options.natomType string = "auto"
        options.nbondType string = "auto"
        options.nangleType string = "auto"
        options.masses string = "auto"
        options.charges string = "none"
        options.dihedrals string = "auto"
        options.ndihedralType string = "auto"
    end

    %%% input
    % atoms - 5 x natom array of molecule IDs, atom types, and positions
    % bonds - 3 x nbond array of bond types and bonded atom IDs
    % angles - 4 xnangle array of angle types and angled atom IDs
    % dihedrals -  4 x ndihedral array of dihedral types and dihedraled atom IDs
    % charges - 1 x natom array of atom charges

    %%% get count for objects
    natom = size(atoms,2);
    nbond = size(bonds,2);
    nangle = size(angles,2);
    nmolecule = max(atoms(1,:));

    %%% interpret input
    if options.natomType == "auto"
        natomType = max(atoms(2,:));
    else
        natomType = str2double(options.natomType);
    end
    if options.nbondType == "auto"
        nbondType = max(bonds(1,:));
    else
        nbondType = str2double(options.nbondType);
    end
    if options.nangleType == "auto"
        nangleType = max(angles(1,:));
    else
        nangleType = str2double(options.nangleType);
    end
    if options.masses == "auto"
        masses = ones(1,natomType);
    else
        masses = options.masses;
    end
    if options.charges == "none"
        charges = zeros(0,natom);
    else
        charges = str2double(options.charges);
        len_charge = floor(log10(max([charges,1])))+1;
    end

    %%% deal with dihedrals
    if options.dihedrals == "auto"
        dihedrals = zeros(5,0);
    else
        dihedrals = str2double(options.dihedrals);
    end
    if options.ndihedralType == "auto"
        ndihedralType = max(dihedrals(1,:));
    else
        ndihedralType = str2double(options.ndihedralType);
    end
    ndihedral = size(dihedrals,2);

    %%% get digit counts
    len_natom = floor(log10(natom))+1;
    len_nbond = floor(log10(nbond))+1;
    len_nangle = floor(log10(nangle))+1;
    len_ndihedral = floor(log10(ndihedral))+1;
    len_nobject = max([len_natom,len_nbond,len_nangle,len_ndihedral]);
    len_nmolecule = floor(log10(nmolecule))+1;
    len_natomType = floor(log10(natomType))+1;
    len_nbondType = floor(log10(nbondType))+1;
    len_nangleType = floor(log10(nangleType))+1;
    len_ndihedralType = floor(log10(ndihedralType))+1;
    len_nobjectType = max([len_natomType,len_nbondType,len_nangleType,len_ndihedralType]);
    len_dbox = floor(log10(dbox/2))+1;

    %%% open file
    f = fopen(geoFile,'w');

    %%% initiate objects
    fprintf(f,"## Number of objects\n");
    fprintf(f,"\t" + ars.fstring(natom,len_nobject,0,"L") + " atoms\n");
    if nbond > 0; fprintf(f,"\t" + ars.fstring(nbond,len_nobject,0,"L") + " bonds\n"); end
    if nangle > 0; fprintf(f,"\t" + ars.fstring(nangle,len_nobject,0,"L") + " angles\n"); end
    if ndihedral > 0; fprintf(f,"\t" + ars.fstring(ndihedral,len_nobject,0,"L") + " dihedrals\n"); end
    fprintf(f,"\n");

    %%% initiate object types
    fprintf(f,"## Number of object types\n");
    fprintf(f,"\t" + ars.fstring(natomType,len_nobjectType,0,"L") + " atom types\n");
    if nbondType > 0; fprintf(f,"\t" + ars.fstring(nbondType,len_nobjectType,0,"L") + " bond types\n"); end
    if nangleType > 0; fprintf(f,"\t" + ars.fstring(nangleType,len_nobjectType,0,"L") + " angle types\n"); end
    if ndihedralType > 0; fprintf(f,"\t" + ars.fstring(ndihedralType,len_nobjectType,0,"L") + " dihedral types\n"); end
    fprintf(f,"\n");

    %%% define simulation box
    fprintf(f,"## Simulation box\n");
    fprintf(f,"\t%0.2f %0.2f xlo xhi\n",-dbox/2,dbox/2);
    fprintf(f,"\t%0.2f %0.2f ylo yhi\n",-dbox/2,dbox/2);
    fprintf(f,"\t%0.2f %0.2f zlo zhi\n\n",-dbox/2,dbox/2);

    %%% define masses
    fprintf(f,"Masses\n\n");
    for i = 1:natomType
        fprintf(f,"\t" + ars.fstring(i,len_natomType,0,"L") + " " + masses(i) + "\n");
    end
    fprintf(f,"\n");

    %%% write atom and their postions
    fprintf(f,"Atoms\n\n");
    for i = 1:natom
        fprintf(f,strcat("\t",...
            ars.fstring(i,len_natom,0,"L"), " ",...
            ars.fstring(atoms(1,i),len_nmolecule,0,"L"), " ",...
            ars.fstring(atoms(2,i),len_natomType,0,"L"), " "));
        if ~isempty(charges)
        fprintf(f,strcat(...
            ars.fstring(charges(i),len_charge+6,4,"R"), " "));
        end
        fprintf(f,strcat(" ",...
            ars.fstring(atoms(3,i),len_dbox+4,2,"R"), " ",...
            ars.fstring(atoms(4,i),len_dbox+4,2,"R"), " ",...
            ars.fstring(atoms(5,i),len_dbox+4,2,"R"), "\n"));
    end
    fprintf(f,"\n");

    %%% write bonds
    if nbond > 0
        fprintf(f,"Bonds\n\n");
        for i = 1:nbond
            fprintf(f,strcat("\t",...
                ars.fstring(i,len_nbond,0,"L"), " ",...
                ars.fstring(bonds(1,i),len_nbondType,0,"L"), "  ",...
                ars.fstring(bonds(2,i),len_natom,0,"L"), " ",...
                ars.fstring(bonds(3,i),len_natom,0,"L"), "\n"));
        end
        fprintf(f,"\n");
    end

    %%% write angles
    if nangle > 0
        fprintf(f,"Angles\n\n");
        for i = 1:nangle
            fprintf(f,strcat("\t",...
                ars.fstring(i,len_nangle,0,"L"), " ",...
                ars.fstring(angles(1,i),len_nangleType,0,"L"), "  ",...
                ars.fstring(angles(2,i),len_natom,0,"L"), " ",...
                ars.fstring(angles(3,i),len_natom,0,"L"), " ",...
                ars.fstring(angles(4,i),len_natom,0,"L"), "\n"));
        end
        fprintf(f,"\n");
    end

    %%% write dihedrals
    if ndihedral > 0
        fprintf(f,"Dihedrals\n\n");
        for i = 1:ndihedral
            fprintf(f,strcat("\t",...
                ars.fstring(i,len_ndihedral,0,"L"), " ",...
                ars.fstring(dihedrals(1,i),len_ndihedralType,0,"L"), "  ",...
                ars.fstring(dihedrals(2,i),len_natom,0,"L"), " ",...
                ars.fstring(dihedrals(3,i),len_natom,0,"L"), " ",...
                ars.fstring(dihedrals(4,i),len_natom,0,"L"), " ",...
                ars.fstring(dihedrals(5,i),len_natom,0,"L"), "\n"));
        end
        fprintf(f,"\n");
    end

    %%% close file
    fclose(f);
end