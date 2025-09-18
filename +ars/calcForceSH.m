%%% simple harmonic force
function F = calcForceSH(r12,r12_eq,k_x)
    F = zeros(3,2);
    F(:,2) = -k_x.*ars.unit_vector(r12)*(norm(r12)-r12_eq);
    F(:,1) = -F(:,2);
end