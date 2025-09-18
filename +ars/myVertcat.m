%%% vertical concatenation of two matrices with same number of columns.
%%% does not throw flag about varaible changing size each loop.
function C = myVertcat(A,B)
    if size(A,2) == size(B,2) || size(A,2) == 0 || size(B,2) == 0
        C = [A;B];
    else
        error("Cannot vertically concatenate matrices with different numbers of columns.")
    end
end