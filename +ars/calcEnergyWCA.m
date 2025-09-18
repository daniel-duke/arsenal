%%% Weeks-Chandler-Anderoson energy
function U = calcEnergyWCA(r12_mag,sigma,epsilon)
    if r12_mag > 0
        U = 4*epsilon*((sigma/r12_mag)^12 - (sigma/r12_mag)^6 + 1/4);
    else
        U = Inf;
    end
end