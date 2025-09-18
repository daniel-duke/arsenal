%%% calculate quaternions for given euler angles
function [chi,eta,xi,zeta] = euler2quat(phi,theta,psi)

%%% Euler angles
% phi    - first rotation (about z-axis)          - azimuthal angle
% theta  - second rotation (about x-axis)         - polar angle
% psi    - third rotation (about rotated z-axis)  - rotor angle

%%% pos quaternion (q = chi + eta*i + xi*j + zeta*k)
chi =  cos(theta/2.0)*cos((psi+phi)/2.0);
eta =  -sin(theta/2.0)*cos((psi-phi)/2.0);
xi  =  sin(theta/2.0)*sin((psi-phi)/2.0);
zeta=  -cos(theta/2.0)*sin((psi+phi)/2.0);

%%% normalization of quaternions
m = sqrt(chi^2 + eta^2 + xi^2 + zeta^2);
chi = chi/m;
eta = eta/m;
xi  = xi/m;
zeta= zeta/m;

end

