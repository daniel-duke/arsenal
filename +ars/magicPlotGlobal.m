%%% set pretty plot properties (applies to everything)
function colors = magicPlotGlobal()
    set(0,'defaultTextInterpreter','latex')
    set(0,'defaultLegendInterpreter','latex')
    set(0,'defaultAxesTickLabelInterpreter','latex')
    set(0,'defaultAxesFontSize',16)
    set(0,'defaultLineLineWidth',1.6);
    set(0,'defaultAxesBox','on');
    colors = lines(7);
end