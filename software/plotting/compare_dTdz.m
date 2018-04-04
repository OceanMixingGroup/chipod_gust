% Compares different temperature gradient estimates.
function [] = compare_dTdz()

    %_____________________include path of processing flies______________________
    addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines

    %____________________set directories______________________
    here    =   pwd;                % mfiles folder
    basedir =   here(1:(end-6));    % substract the mfile folder
    savedir =   [basedir 'proc/'];  % directory directory to save data
    unit    = chi_get_unit_name(basedir); % get unit name
    rawdir       = [basedir filesep 'raw' filesep]; % raw files location

    %_____________________set time limits______________________
    % get time limits from whoAmI;
    [TL] = whoAmI_timeLimits(basedir);
    time_lim = TL.master;

    if exist([basedir filesep 'input' filesep 'dTdz_m.mat'], 'file')
        load ../input/dTdz_m.mat;
    end

    if exist([basedir filesep 'input' filesep 'dTdz_i.mat'], 'file')
        load ../input/dTdz_i.mat;
    end

    if exist([basedir filesep 'input' filesep 'dTdz_w.mat'], 'file')
        load ../input/dTdz_w.mat;
    end

    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
                 'Papersize',[30 20],'PaperPosition',[0 0 30 20]);

    [ax, ~] = create_axes(fig, 3, 1, 0);

    tl = time_lim;

    a = 1;
    if isfield(Tz_i,'T')
        hdl = plot(ax(a), Tz_i.time, Tz_i.T, 'Linewidth', 1);
        xlim(ax(a), tl);

        t = text_corner(ax(a), ['Temperature [^\circ C]'], 1);
        t.Color = hdl.Color;
    else
        hdl1 = plot(ax(a), Tz_i.time, Tz_i.T1, 'Linewidth', 1);
        hdl2 = plot(ax(a), Tz_i.time, Tz_i.T2, 'Linewidth', 1);

        t = text_corner(ax(a), ['T1 [^\circ C]'], 1);
        t.Color = hdl1.Color;
        t = text_corner(ax(a), ['T2 [^\circ C]'], 3);
        t.Color = hdl2.Color;

        % plot(ax(a), Tz_i.time, Tz_i.T12, 'Linewidth', 1);
        xlim(ax(a), tl);
    end

    a = 2;
    if isfield(Tz_i,'T')
        po10 = floor(log10(max(abs(Tz_i.Tz))));
        hdl1 = plot(ax(a), Tz_i.time, Tz_i.Tz/10^po10, 'Linewidth', 1);
        t = text_corner(ax(a), ['Tz_i [^\circ C]'], 1);
        t.Color = hdl1.Color;

        if exist('Tz_w', 'var')
            ndt = round(diff(Tz_i.time(1:2))/diff(Tz_w.time(1:2)));
            hdl2 = plot(ax(a), ...
                        moving_average(Tz_w.time, ndt, ndt), ...
                        moving_average(Tz_w.Tz, ndt, ndt)/10^po10, 'linewidth', 1);

            t = text_corner(ax(a), ['Tz_{sorted} [^\circ C]'], 3);
            t.Color = hdl2.Color;
        end
    else
        po10 = floor(log10(max(abs(Tz_i.Tz1))));
        hdl1 = plot(ax(a), Tz_i.time, Tz_i.Tz1/10^po10, 'Linewidth', 1);
        hdl2 = plot(ax(a), Tz_i.time, Tz_i.Tz2/10^po10, 'Linewidth', 1);
        % plot(ax(a), Tz_i.time, Tz_i.Tz12/10^po10, 'Linewidth', 1);

        t = text_corner(ax(a), ['Tz^1_{internal} [10^{' num2str(po10) '}K/m]'], 1);
        t.Color = hdl1.Color;

        t = text_corner(ax(a), ['Tz^2_{internal} [10^{' num2str(po10) '}K/m]'], 3);
        t.Color = hdl2.Color;

        if exist('Tz_w', 'var')
            ndt = round(diff(Tz_i.time(1:2))/diff(Tz_w.time(1:2)));
            hdl1 = plot(ax(a), ...
                        moving_average(Tz_w.time, ndt, ndt), ...
                        moving_average(Tz_w.Tz1, ndt, ndt)/10^po10, 'linewidth', 1);
            hdl2 = plot(ax(a), ...
                        moving_average(Tz_w.time, ndt, ndt), ...
                        moving_average(Tz_w.Tz2, ndt, ndt)/10^po10, 'linewidth', 1);

            t = text_corner(ax(a), ['Tz^1_{sorted} [10^{' num2str(po10) '}K/m]'], 7);
            t.Color = hdl1.Color;

            t = text_corner(ax(a), ['Tz^2_{sorted} [10^{' num2str(po10) '}K/m]'], 9);
            t.Color = hdl2.Color;
        end

        if exist('Tz_m', 'var')
            hdl = plot(ax(a), Tz_m.time, Tz_m.Tz/10^po10, 'k', 'linewidth', 1);

            t = text_corner(ax(a), ['Tz^2_{moor} [10^{' num2str(po10) '}K/m]'], 8);
            t.Color = hdl.Color;
        end
    end
    plot(ax(a), tl, [0 0],':k', 'Linewidth', 1);
    xlim(ax(a), tl);

    a = 3;
    if isfield(Tz_i,'T')
        po10 = floor(log10(max(abs(Tz_i.N2))));
        plot(ax(a), Tz_i.time, Tz_i.N2/10^po10, 'Linewidth', 1);
    else
        po10 = floor(log10(max(abs(Tz_i.N2_1))));
        plot(ax(a), Tz_i.time, Tz_i.N2_1/10^po10, 'Linewidth', 1);
        plot(ax(a), Tz_i.time, Tz_i.N2_2/10^po10, 'Linewidth', 1);
        plot(ax(a), Tz_i.time, Tz_i.N2_12/10^po10, 'Linewidth', 1);
    end
    plot(ax(a), tl, [0 0],':k', 'Linewidth', 1);
    xlim(ax(a), tl);
    t = text_corner(ax(a), ['N^2 [10^{' num2str(po10) '} s^{-2}]'], 1);

    datetick(ax(a), 'keeplimits');

    t = text_corner(ax(1), ['T_z of unit ' unit], -2);

    linkaxes(ax, 'x')
    xlim(ax(1), time_lim);

    print(gcf,'../pics/Compare_dTdz.png','-dpng','-r200','-painters');
    savefig(fig, '../pics/Compare_dTdz.fig');

end