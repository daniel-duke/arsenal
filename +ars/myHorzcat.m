%%% horizontal concatenation of two matrices with same number of rows.
%%% does not throw flag about varaible changing size each loop.
function C = myHorzcat(A,B)
    if size(A,1) == size(B,1) || size(A,1) == 0 || size(B,1) == 0
        C = [A,B];
    else
        error("Cannot horizontally concatenate matrices with different numbers of rows.")
    end
end