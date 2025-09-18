%%% bending energy
function U = calcEnergyBend(r12,r23,theta_eq,k_theta)
    C32 = sum(r23.*r12);
    cos = C32/(norm(r12)*norm(r23));
    theta = acos(cos);
    U = 1/2*k_theta*(theta-theta_eq)^2;
end