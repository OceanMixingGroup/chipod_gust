function [ax] = DebugRawData(T, t0, t1, Turb, nantimes)

    tind = find_approx(T.time, t0, 1):find_approx(T.time, t1, 1);

    if exist('nantimes', 'var')
        for sensor=1:3
            if ~isempty(nantimes{sensor})
                for tt=1:size(nantimes{sensor}, 1)
                    i0 = find_approx(T.time(tind), nantimes{sensor}(tt, 1));
                    i1 = find_approx(T.time(tind), nantimes{sensor}(tt, 2));

                    nanind{sensor}(tt,:) = [i0, i1];
                end
            else
                nanind{sensor} = [];
            end
        end
    end

    ax(1) = subplot(221);
    PlotVar(T.time(tind), T.T1Pt(tind) - nanmean(T.T1Pt(tind)) + 0.1, nanind{1});
    hold on;
    ax(1).ColorOrderIndex = 2;
    PlotVar(T.time(tind), T.T2Pt(tind) - nanmean(T.T2Pt(tind)), nanind{2});
    ylim([-0.1 0.2])
    datetick('keeplimits');
    legend('T1P', 'T2P')
    ylabel('TP')

    ax(2) = subplot(222);
    semilogy(moving_average(T.time(tind), 60, 30), moving_var(T.T1Pt(tind), 60, 30));
    hold on;
    semilogy(moving_average(T.time(tind), 60, 30), moving_var(T.T2Pt(tind), 60, 30));
    ylim([1e-7, 1e-3])
    datetick('keeplimits');
    legend('T1P', 'T2P')
    ylabel('var(TP)')

    fnames = fieldnames(Turb);
    for ff = 1:length(fnames)
        if ~strcmpi(fnames{ff}(1:3), 'chi') | strcmpi(fnames{ff}(end-1:end), 'ic')
            continue;
        end

        chi = Turb.(fnames{ff});

        Turbtind = find_approx(chi.time, t0, 1):find_approx(chi.time, t1, 1);

        ax(4) = subplot(224);
        semilogy(chi.time(Turbtind), chi.chi(Turbtind), 'displayname', fnames{ff}(5:end))
        hold on;
    end

    legend(ax(4))

    ax(3) = subplot(223);
    [axx, h1, h2] = plotyy(T.time(tind), T.W(tind), ...
                           chi.time(Turbtind), chi.dTdz(Turbtind));
    h1.delete;
    hold(axx(2), 'on')
    hold(ax(3), 'on')
    axx(1).ColorOrderIndex = 1;
    PlotVar(T.time(tind), T.W(tind), nanind{3}, axx(1))
    h2.LineWidth = 2;
    h2.LineStyle = '-';
    datetick('keeplimits');
    ylabel(axx(1), 'W')
    ylabel(axx(2), 'dT/dz')
    ylim(axx(2), [-1 1]*7e-3);
    axx(2).YTick = [-7, -5, -3, -1, 1, 3, 5, 7] * 1e-3;
    plot(axx(2), xlim(axx(2)), [1 1]*1e-3, 'k--', ...
         xlim(axx(2)), [1 1]*-1e-3, 'k--', ...
         xlim(axx(2)), [0 0], 'k--');

    % 0-crossings of dT/dz
    sgn = sign(chi.dTdz(Turbtind));
    tind0cross = Turbtind(abs(diff(sgn)) > 0);
    Mark0Cross(ax, chi.time(tind0cross));

    axes(ax(1)); title('minor ticks = 0-crossings of dT/dz');

    axes(ax(4));
    ylim([1e-10, 1e-3]);
    datetick('keeplimits')
    ylabel('\chi')

    linkaxes([ax, axx(2)], 'x');
    xlim([t0, t1])

    for aa = 1:4
        set(ax(aa), 'Units', 'normalized')
        ax(aa).Position(1) = ax(aa).Position(1) - 0.05;
        ax(aa).Position(2) = ax(aa).Position(2) - 0.05;
        ax(aa).Position(3) = ax(aa).Position(3) + 0.08;
        ax(aa).Position(4) = ax(aa).Position(4) + 0.08;
    end
    ax(2).Position(1) = ax(2).Position(1) + 0.05;
    ax(4).Position(1) = ax(4).Position(1) + 0.05;

end

function [] = PlotVar(time, var, naninds, ax)

    if ~exist('ax', 'var'), ax = gca(); end

    plot(ax, time, var);
    hold on;
    for tt=1:size(naninds, 1)
        nantind = naninds(tt,1):naninds(tt,2);
        plot(ax, time(nantind), var(nantind), '-',...
             'color', [1 1 1]*0.6, 'HandleVisibility', 'off');
    end
end

function [] = Mark0Cross(ax, tind)
    for aa = 1:length(ax)
        ax(aa).XAxis.TickDirection= 'in';
        ax(aa).XMinorTick = 'on';
        ax(aa).XAxis.MinorTickValues = tind;
        ax(aa).XAxis.TickLength(1) = 0.08;
    end

    % for tt=1:length(tind0cross)
    %     for ii=[1, 4]
    %         plot(ax(ii), [1, 1]*chi.time(tind0cross(tt)), ax(ii).YLim, ...
    %              '-', 'color', [1 1 1]*0.5, 'HandleVisibility', 'off');
    %     end
    % end

end