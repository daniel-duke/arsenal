%%% correlation between two vectors
function v = calcUnitDot(v1,v2)
    v1_u = ars.unitVector(v1);
    v2_u = ars.unitVector(v2);
    v = dot(v1_u,v2_u);
end