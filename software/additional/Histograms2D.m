function Histograms2D(chi, ID)

    dt = (chi.time(2)-chi.time(1))*86400;

    CreateFigure;
    ax(1) = subplot(321);
    [counts, xbins, ybins] = hc(chi.dTdz, log10(chi.chi));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('dT/dz'); ylabel('log_{10}\chi')
    title([ID ' | ' num2str(dt) 's estimates'])

    ax(3) = subplot(323);
    [counts, xbins, ybins] = hc(chi.dTdz, log10(chi.eps));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('dT/dz'); ylabel('log_{10}\epsilon')

    ax(2) = subplot(322);
    [counts, xbins, ybins] = hc(chi.spd, log10(chi.chi));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('spd'); ylabel('log_{10}\chi')

    ax(4) = subplot(324);
    [counts, xbins, ybins] = hc(chi.spd, log10(chi.eps));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('spd'); ylabel('log_{10}\epsilon')

    ax(5) = subplot(325);
    [counts, xbins, ybins] = hc(chi.dTdz, log10(chi.Kt));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('dTdz'); ylabel('log_{10}K_T')

    ax(6) = subplot(326);
    [counts, xbins, ybins] = hc(chi.spd, log10(chi.Kt));
    pcolor(xbins, ybins, counts); shading flat;
    xlabel('spd'); ylabel('log_{10}K_T')

    linkaxes(ax(1:2), 'y');
    linkaxes(ax(3:4), 'y');

    linkaxes(ax([1, 3, 5]), 'x');
    linkaxes(ax([2, 4, 6]), 'x');

    colormap(flipud(gray))
end

function [counts, xbins, ybins] = hc(x, y)

    mask = ~isnan(x) & ~isnan(y);
    [counts, xe, ye] = histcounts2(x(mask), y(mask), ...
                                   'numbins', round(sqrt([numel(x) numel(y)])));

    counts = counts'./max(counts(:));

    % try histogram equalization
    try
        counts = histeq(counts, 200);
    catch ME
    end

    xbins = (xe(1:end-1) + xe(2:end))/2;
    ybins = (ye(1:end-1) + ye(2:end))/2;
end