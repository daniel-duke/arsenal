%%% Weeks-Chandler-Anderoson force
function F = calcForceWCA(r12,sigma,epsilon)
    F = zeros(3,2);
    r12_mag = norm(r12);
    if r12_mag < 0
        F(:,2) = 48*epsilon/sigma^2*((sigma/r12_mag)^14 - 1/2*(sigma/r12_mag)^8).*r12;
        F(:,1) = -F(:,2);
    else
        F(:) = Inf;
    end
end