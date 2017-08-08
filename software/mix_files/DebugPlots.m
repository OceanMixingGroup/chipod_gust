function [] = DebugPlots(hfig, t0, t1, chi, name, ww)

    if ~exist('ww', 'var'), ww = 1; end
    if isempty(hfig), hfig = gcf(); end

    i0 = find_approx(chi.time, t0, 1);
    i1 = find_approx(chi.time, t1, 1);
    tind = i0:i1;

    figure(hfig);
    ax(1) = subplot(511);
    semilogy(chi.time(tind), chi.N2(tind)./chi.dTdz(tind).^2, 'displayname', name)
    ylabel('N^2/T_z^2')
    Common()

    ax(2) = subplot(512);
    plot(chi.time(tind), chi.dTdz(tind), 'displayname', name)
    hold on;
    plot(xlim, [0, 0], 'k--');
    ylabel('dT/dz')
    Common()

    ax(3) = subplot(513);
    try
        semilogy(chi.time(tind), chi.eps(tind), 'displayname', name)
    catch ME
        semilogy(chi.time(tind), chi.eps1(tind), 'displayname', name)
    end
    ylabel('\epsilon')
    Common()

    ax(4) = subplot(514);
    try
        semilogy(chi.time(tind), chi.Kt(tind), 'displayname', name)
    catch ME
        semilogy(chi.time(tind), chi.Kt1(tind), 'displayname', name)
    end
    ylabel('K_t')
    Common()

    ax(5) = subplot(515);
    try
        plot(moving_average(chi.time(tind), ww, ww), ...
             moving_average(chi.Jq(tind), ww, ww), 'displayname', name)
    catch ME
        plot(moving_average(chi.time(tind), ww, ww), ...
             moving_average(-chi.Jq1(tind), ww, ww), 'displayname', name)
    end
    ylabel('J_q^t')
    Common()
    legend('-DynamicLegend');

    linkaxes(ax, 'x')
    xlim([t0, t1])
    datetick('x', 'mmm-dd HH:MM', 'keeplimits')
end

function Common()
    hold on
    datetick('x', 'mm-dd HH:MM', 'keeplimits')
end