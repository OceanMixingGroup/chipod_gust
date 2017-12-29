function [hfig] = CreateFigure(vis)
    if nargin <1
        vis = 'on';
    end

    hfig = figure('Color',[1 1 1],'visible', vis, ...
                  'Position', [100 100 1400 900]);

    set(hfig, 'DefaultAxesFontSize', 18);
    set(hfig, 'renderer', 'opengl');

    set(hfig,'DefaultAxesBox','on')
    set(hfig,'DefaultAxesTickDir','out')

    set(hfig, 'DefaultHistogramEdgeColor', 'none');

end