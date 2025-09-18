%%% turn vector into unit vector
function v = unitVector(v)
    v = v./norm(v);
end