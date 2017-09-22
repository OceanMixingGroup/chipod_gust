function Histograms(chi, hfig, normstr, ID)

    nbins = 300;
    chibins = linspace(-10, -3, nbins);
    epsbins = linspace(-12, -2, nbins);
    ktbins  = linspace(-8, 0, nbins);
    jqbins  = linspace(-3, 4, nbins/2);

    figure(hfig)
    hfig.Position(3) = 2000;
    hfig.Position(4) = 1065;
    set(hfig, 'DefaultLegendBox', 'off')

    hax(1) = subplot(221);
    set(gca, 'color', 'none')
    myhist(chi.chi, chibins, normstr, ID)
    hold on;
    xlabel('log_{10} \chi')
    xlim([-12, -4])
    ylabel(normstr)

    hax(2) = subplot(222);
    set(gca, 'color', 'none')
    myhist(chi.eps, epsbins, normstr, ID)
    xlim([-12, -2])
    xlabel('log_{10} \epsilon')
    ylabel(normstr)
    hold on;

    hax(3) = subplot(223);
    set(gca, 'color', 'none')
    myhist(chi.Kt, ktbins, normstr, ID)
    xlim([-8, 2])
    hold on;
    xlabel('log_{10} K_T')
    ylabel(normstr)

    hax(4) = subplot(224);
    set(gca, 'color', 'none')
    myhist(abs(chi.Jq), jqbins, normstr, ID)
    xlim([-5, 4.5])
    hold on;
    xlabel('log_{10} |J_q|')
    ylabel(normstr)
end

function myhist(var, bins, normstr, ID)

    avg = nanmean(var);
    med = nanmedian(var);

    ms = 8; % markersize

    str = sprintf('%s\n%s', ID, ...
                  ['\mu=' num2str(avg, '%.1e') '|mdn=' num2str(med, '%.1e')]);

    color = choose_color(ID,'color');
    lw = choose_color(ID,'width');

    hh = histogram(log10(var), 'BinEdges', bins, 'normalization', normstr, ...
                   'displayname', str, 'displaystyle', 'stairs', ...
                   'LineWidth', lw, 'EdgeColor', color);

    if strcmpi(normstr, 'pdf')
        ylim([0 max([0.5, ylim])])
    end
    hold on;

    ylims = ylim;
    plot(log10(avg), ylims(2)*0.9, 'v', 'color', hh.EdgeColor, ...
         'markersize', ms, 'handlevisibility', 'off');

    plot(log10(med), ylims(2)*0.9 , '+', 'color', hh.EdgeColor, ...
         'markersize', ms, 'handlevisibility', 'off');

    legend('-dynamiclegend', 'Location', 'northeastoutside');
end
