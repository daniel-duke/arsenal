%%% three perpendicular unit vectors
function R = randBasis()
    z_basis = ars.randUnitVec();
    y_basis = ars.unitVector(cross(z_basis,ars.randUnitVec()));
    x_basis = cross(y_basis,z_basis);
    R = [x_basis,y_basis,z_basis];
end