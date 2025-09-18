%%% simple harmonic energy
function U = calcEnergySH(r12,r12_eq,k_x)
    U = 1/2*k_x*(norm(r12)-r12_eq)^2;
end