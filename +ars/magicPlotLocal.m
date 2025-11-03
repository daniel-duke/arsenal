%%% set pretty plot properties (applies to given figure)
function colors = magicPlotLocal(fig)
    arguments
        fig matlab.ui.Figure = gcf
    end
    set(fig,'defaultTextInterpreter','latex')
    set(fig,'defaultLegendInterpreter','latex')
    set(fig,'defaultAxesTickLabelInterpreter','latex')
    set(fig,'defaultAxesFontSize',16)
    set(fig,'defaultLineLineWidth',2);
    set(fig,'defaultConstantLineLineWidth',2);
    set(fig,'defaultAxesBox','on');
    colors = lines(7);
end