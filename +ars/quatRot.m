%%% pos quaternion rotation matrix
function A = quatRot(chi,eta,xi,zeta)
    A(1,1) = -zeta^2 + eta^2 - xi^2 + chi^2;
    A(1,2) = 2*(xi*eta - zeta*chi);
    A(1,3) = 2*(eta*zeta + xi*chi);
    
    A(2,1) = 2*(xi*eta + zeta*chi);
    A(2,2) = -zeta^2 - eta^2 + xi^2 + chi^2;
    A(2,3) = 2*(xi*zeta - eta*chi);
    
    A(3,1) = 2*(eta*zeta - xi*chi);
    A(3,2) = 2*(xi*zeta + eta*chi);
    A(3,3) = zeta^2 - eta^2 - xi^2 + chi^2;
end