%%% calculate center of mass, using method from Bai and Breen 2008
function com = calcCOM(r,dbox)
    xi_bar = mean(cos(2*pi.*(r./dbox + 1/2)),2);
    zeta_bar = mean(sin(2*pi.*(r./dbox + 1/2)),2);
    theta_bar = atan2(-zeta_bar,-xi_bar) + pi;
    r_ref = dbox.*(theta_bar./(2*pi) - 1/2);
    com = r_ref + mean(ars.applyPBC(r - r_ref, dbox),2);
end