%%% set pretty colorbar properties
function c = magicColorbar(ax)
    arguments
        ax matlab.graphics.axis.Axes = gca
    end
    c = colorbar;
    c.Label.Interpreter = ax.TickLabelInterpreter;
    c.TickLabelInterpreter = ax.TickLabelInterpreter;
    c.FontSize = ax.FontSize;
end