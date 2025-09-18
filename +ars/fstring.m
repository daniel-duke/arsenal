%%% write number in string form with given parameters
function fstr = fstring(x,space,precision,alignment,fill)
    arguments
        x double
        space double
        precision double = 0
        alignment string = "R"
        fill string = "space"
    end
    if space == 0
        fstr = sprintf(sprintf("%%0.%df",precision),x);
    elseif alignment == "L"
        fstr = sprintf(sprintf("%%-%d.%df",space,precision),x);
    elseif alignment == "R"
        if fill == "space"
            fstr = sprintf(sprintf("%%%d.%df",space,precision),x);
        elseif fill == "zero"
            fstr = sprintf(sprintf("%%0%d.%df",space,precision),x);
        end
    end
end