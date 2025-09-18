%%% test position against list of other positions for overlap
function overlap = checkOverlap(r,r_other,sigma,dbox)
    overlap = false;
    for i = 1:size(r_other,2)
        if norm(ars.applyPBC(r_other(:,i) - r, dbox)) < sigma
            overlap = true;
            return
        end
    end
end