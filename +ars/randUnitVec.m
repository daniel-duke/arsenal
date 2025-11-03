%%% unit vector with random orientation and given number of dimensions
function x = randUnitVec(ndim)
    arguments
        ndim double = 3
    end
    x = ars.unitVector(randn(ndim,1));
end