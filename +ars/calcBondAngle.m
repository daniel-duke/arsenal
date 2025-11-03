%%% angle between two bond vectors, range [0,180]
function theta = calcBondAngle(v1,v2)
    theta = 180 - acosd(ars.calcUnitDot(v1,v2));
end