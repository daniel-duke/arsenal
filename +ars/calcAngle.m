%%% angle between two vectors, range [0,180]
function theta = calcAngle(v1,v2)
    theta = acosd(ars.calcUnitDot(v1,v2));
end