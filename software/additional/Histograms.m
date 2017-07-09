function Histograms(chi, hfig, normstr, legstr)

    subplot(221);
    set(gca, 'color', 'none')
    myhist(chi.chi, normstr, legstr)
    hold on;
    xlabel('log_{10} \chi')
    ylabel(normstr)

    subplot(222);
    set(gca, 'color', 'none')
    myhist(chi.eps, normstr, legstr)
    xlabel('log_{10} \epsilon')
    ylabel(normstr)
    hold on;

    subplot(223)
    set(gca, 'color', 'none')
    myhist(chi.Kt, normstr, legstr)
    hold on;
    xlabel('log_{10} K_T')
    ylabel(normstr)

    subplot(224)
    set(gca, 'color', 'none')
    myhist(abs(chi.Jq), normstr, legstr)
    hold on;
    xlabel('log_{10} |J_q|')
    ylabel(normstr)
end

function myhist(var, normstr, legstr)
    histogram(log10(var), 'normalization', normstr, ...
              'binmethod', 'sqrt', 'displayname', legstr, ...
              'edgecolor', 'none', 'facealpha', 0.4, ...
              'displaystyle', 'stairs')
end
