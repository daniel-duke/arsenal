%%% get random unit vector with given number of dimensions
function x = unitVecRand(ndim)
    arguments
        ndim double = 3
    end
    x = ars.unitVector(randn(ndim,1));
end