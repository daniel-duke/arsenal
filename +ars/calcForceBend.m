%%% bending force
function F = calcForceBend(r12,r23,theta_eq,k_theta)
    F = zeros(3,3);
    C22 = sum(r12.*r12);
    C33 = sum(r23.*r23);
    C32 = sum(r23.*r12);
    C3322SQ = sqrt(C33*C22);
    cos = C32/C3322SQ;
    theta = acos(cos);
    theta_eff = theta-theta_eq;
    if theta_eff ~= 0
        dCos3 = -(C32/C33.*r23 - r12)/C3322SQ;
        dCos2 =  (C32/C33.*r23 - C32/C22.*r12 + r23 - r12)/C3322SQ;
        dCos1 =  (C32/C22.*r12 - r23)/C3322SQ;
        dThdCos = -1/sqrt(1-cos^2);
        F(:,3) = -k_theta*dThdCos*dCos3*theta_eff;
        F(:,2) = -k_theta*dThdCos*dCos2*theta_eff;
        F(:,1) = -k_theta*dThdCos*dCos1*theta_eff;
    end
end