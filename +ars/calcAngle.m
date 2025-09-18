%%% angle between two vectors
function v = calcAngle(v1,v2)
    v = acosd(ars.calcUnitDot(v1,v2));
end