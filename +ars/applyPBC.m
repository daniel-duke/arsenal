%%% apply periodic boundary condition to a position vector
function rPBC = applyPBC(r,dbox)
    rPBC = r - dbox.*round(r./dbox);
end