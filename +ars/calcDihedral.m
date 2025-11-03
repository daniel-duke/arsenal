%%% dihedral angle between three bond vectors, range (-180,180]
function phi = calcDihedral(v1,v2,v3)
    v2 = ars.unitVector(v2);
    n1 = cross(v1,v2);
    n2 = cross(v2,v3);
    x = dot(n1,n2);
    y = dot(v2,cross(n1,n2));
    phi = atan2d(y,x);
end