% [wda] = winters_dasaro_avg(t0, t1, vdisp, chi, T, Tp, dt, plotflag)
% Treats the chipod as a profiler and applies the method of Winters & D'Asaro to estimate heat flux.
% Inputs:
%        t0, t1 - time indices for accelerometer data (usually 50Hz)
%        vdisp, chi, T, Tp - structures passed in from do_wda_estimate
%        plotflag - if 1, make informative plot showing resorted profile.
%        dt - pflag.master.wda_dt (time chunk over which to sort and average)
% Outputs:
%        wda - structure with fields
%           tstart, tstop - start,end of time chunk
%

function [wda] = winters_dasaro_avg(t0, t1, vdisp, chi, T, Tp, dt, plotflag)

    optional = 0;

    MIN_DIS_Z = 0.2; % minimum length (in metres) of a single up- or down-cast
    nquantiles = 10; round(2/60 * dt); % effectively number of bins

    chi_is_empty = ~isfield(chi, 'chi');

    % determine "profile" boundaries
    [~, locs1] = findpeaks(vdisp.dis_z(t0:t1));
    [~, locs2] = findpeaks(-vdisp.dis_z(t0:t1));
    locs = sort([t0; locs1+t0-1; locs2+t0-1; t1]);

    chit0 = find_approx(chi.time, vdisp.time(t0));
    chit1 = find_approx(chi.time, vdisp.time(t1));

    zfull = -(vdisp.dis_z(t0:t1)-nanmean(vdisp.dis_z(t0:t1))); % t0 element can be nan?
    Tfull = T.Tenh(t0:t1);
    Tchi = chi.T(chit0:chit1);

    tstart = chi.time(chit0);
    tstop = chi.time(chit1);

    % initialized returned structure
    wda.nbins = nquantiles+2;
    wda.Tbins = (nan(nquantiles+2, 1));
    wda.tstart = tstart;
    wda.tstop = tstop;
    wda.dTdz_bins = (nan(nquantiles+2, 1));
    wda.dz = (nan(nquantiles+2, 1));
    wda.no_min_dz = 0;

    if all(isnan(T.Tenh)), return; end

    % fundamental assumption used below
    assert(isequal(Tp.time, vdisp.time));

    % loop over valid profiles and save 1 Hz temp observations in those profiles.
    % these 1Hz observations are used to identify isothermal surfaces Tbins.
    il = 1; Tprof = [];
    for ll=1:length(locs)-1
        l0 = locs(ll); l1 = locs(ll+1);
        if abs(vdisp.dis_z(l0) - vdisp.dis_z(l1)) > MIN_DIS_Z
            lstart(il) = locs(ll);
            lstop(il) = locs(ll+1);
            il = il+1;

            ct0 = find_approx(chi.time(chit0:chit1), vdisp.time(l0)) + chit0-1;
            ct1 = find_approx(chi.time(chit0:chit1), vdisp.time(l1)) + chit0-1;

            Tsubset = chi.T(ct0:ct1);

            % only look at those temperature measurements with valid chi estimates
            % chisubset = chi.chi(ct0:ct1);
            % chisubset(chi.spec_area(ct0:ct1) < 2 * chi.spec_floor * chi.nfft) = 0;
            % Tsubset = Tsubset(~isnan(chisubset));

            % save 1sec avg temp in valid profiles.
            % This is *only* used to make bins
            Tprof = [Tprof Tsubset];
        end
    end

    % Lets make some T bins by splitting into quantiles
    % so that each bin has equal number of observations
    Tbins = generate_wda_bins(unique(round(Tprof, 5)), nquantiles);
    dT = diff(Tbins);

    % if no long enough profiles are found
    if isempty(Tprof), wda.no_min_dz = 1; return; end

    % accelerometer has crapped out and thinks we've covered >10m in a minute
    if max(zfull) - min(zfull) > 10, return; end

    zthorpe = linspace(min(zfull), max(zfull), 1000);
    if optional
        Tj = linspace(min(Tfull), max(Tfull), 1000);
        zTj = nan(length(Tj), length(locs));
    end
    dzmat = nan(length(dT), length(locs));

    % set up figure + plot timeseries
    if plotflag

        CreateFigure;
        gray = [1, 1, 1] * 0.75;

        % temp profiles
        htemp = subplot(4,4,[4, 8, 12]);
        hold on; xlabel('T [C]'); title('Unsorted profiles');
        htemp.YTickLabels = [];

        % sorted temp profiles
        hsort = subplot(4,4,[3, 7, 11]); hold on; ylabel('z [m]')
        % hsort.XGrid = 'on';
        hsort.XAxis.TickDirection = 'out';
        hsort.XTick = Tbins;
        hsort.XTickLabels = [];
        hsort.XAxis.TickLabelRotation = 30; hsort.GridAlpha = 1;
        title(hsort, 'Individual sorted profiles')

        % displacement
        hdisp = subplot(4,4,[1, 2]);
        hzfull = plot(hdisp, vdisp.time(t0:t1), -(vdisp.dis_z(t0:t1)), ...
                      'color', gray, 'linewidth', 2);
        hold on;
        xlabel('time'); ylabel('Displacement z [m]')
        hdisp.Color = 'none';
        % enhanced temp time series

        % htemp2 = axes('Position', hdisp.Position, 'Color', 'none'); hold on;
        htemp2 = subplot(4,4,[5, 6]);
        tpt0 = find_approx(Tp.time, vdisp.time(t0));
        tpt1 = find_approx(Tp.time, vdisp.time(t1));
        hfullT = plot(htemp2, T.time_enh(tpt0:tpt1), T.Tenh(tpt0:tpt1), ...
                      'color', gray);
        ylim(htemp2, robust_lim(T.Tenh(tpt0:tpt1)))
        hold on; xlabel('time'); ylabel('T [C]')

        % Tp time series
        htp = subplot(4,4,[9, 10]); hold on;  xlabel('time'); ylabel('dT/dt [C/s]')
        htpfull = semilogy(htp, Tp.time(tpt0:tpt1), Tp.tp(tpt0:tpt1), ...
                           'color', gray);
        htp.Clipping = 'off';
        % ylim(htp, robust_lim(Tp.tp(tpt0:tpt1)))

        if ~chi_is_empty
            hchi = subplot(4,4,[13, 14]);
            hcfull = plot(hchi, chi.time(chit0:chit1), ...
                          chi.chi(chit0:chit1), 'HandleVisibility', 'off', ...
                          'color', gray, 'linewidth', 2);
            hchi.YScale = 'log';
            hold on;
            xlabel(['time since ' datestr(chi.time(chit0)) ' [MM:ss]']);
            ylabel(['\chi [C^2/s]'])

            hchi_scatter = subplot(4,4,15);
            hscfull = scatter(hchi_scatter, ...
                              Tchi, chi.chi(chit0:chit1), ...
                              32, repmat(gray, [length(Tchi), 1]), ...
                              'filled', 'HandleVisibility', 'off');

            hold on; xlabel('T [C]')
            hchi_scatter.XGrid = 'on';
            hchi_scatter.YScale = 'log';
            hchi_scatter.XAxis.TickDirection = 'out';
            hchi_scatter.XAxis.TickLabelRotation = 90;
            hchi_scatter.XTick = Tbins;
            hchi_scatter.GridAlpha = 1;
            fs = hchi_scatter.FontSize;
            hchi_scatter.FontSize = fs * 0.8;

            linkaxes([hsort, hchi_scatter], 'x')
            linkaxes([hchi_scatter, hchi], 'y')

            hchi.YLim = robust_lim(chi.chi(chit0:chit1));
        end

        % hcloud = plot(hsort, Tfull, zfull, '.', ...
        %               'color', [1 1 1]*0.8, 'HandleVisibility', 'off');
        % uistack(hcloud, 'bottom');

        % noise floor diagnostics
        % if ~chi_is_empty
        %     semilogy(hchi, chi.time(chit0:chit1), chi.spec_area(chit0:chit1), ...
        %              'color', 'k', 'linewidth', 2, 'displayname', ['spec_area'])
        %     legend(hchi, '-dynamiclegend')
        %     plot(hchi, hchi.XLim, [1, 1]*chi.spec_floor * chi.nfft, 'k-', ...
        %          'displayname', 'noise floor');
        %     plot(hchi, hchi.XLim, 2*[1, 1]*chi.spec_floor * chi.nfft, 'k--', ...
        %          'displayname', '2x noise floor');
        % end

        hleg = legend(htemp, '-dynamiclegend');
        hleg.Position = [0.734 0.132 0.176 0.092];
        hleg.Title.String = 'Various T_z estimates';
        hleg.Box = 'off';

        % hleg.Position = [0.7479    0.7733    0.1650    0.0800];

        for aa=[hdisp, htemp2, htp]
            datetick(aa, 'x', 'MM:SS', 'keeplimits');
        end

        if ~chi_is_empty
            linkaxes([hchi hdisp htp htemp2], 'x');
            datetick(hchi, 'x', 'MM:SS', 'keeplimits');
            axx = [hdisp, htemp2, htp, hchi, hsort, hchi_scatter, htemp];
        else
            axx = [hdisp, htemp2, htp, hsort, htemp];
            linkaxes([hdisp htp htemp2], 'x');
        end

        x0 = 0.05; y0=0.9;
        labels = 'abcdefgh';
        for xxx=1:length(axx)
            text(axx(xxx), x0, y0, ['(' labels(xxx) ')'], 'units', 'normalized')
        end

        % set reasonable limits
        hdisp.XLim = chi.time([chit0, chit1]);
        hsort.XLim = htemp2.YLim;

        yloc = NaN;
    end

    for ll=1:length(lstart) % loop over "profiles"
        l0 = lstart(ll); l1 = lstop(ll);

        if abs(vdisp.dis_z(l0) - vdisp.dis_z(l1)) < MIN_DIS_Z, continue; end

        zvec = -vdisp.dis_z(l0:l1);
        Tvec = T.Tenh(l0:l1); %T.T(l0:l1);

        if all(isnan(Tvec)), continue; end

        % Tp glitch; ignore this profile
        if any(Tp.tp(l0:l1) > 15), continue; end

        % interpolate to uniform depths before sorting
        Tinterp = interp1(zvec, Tvec, zthorpe);

        % At least 5 points in profile
        if sum((~isnan(Tinterp))) < 5, continue; end

        % Tsort is on depth grid zthorpe
        Tsort = thorpeSort(Tinterp);

        % remove the tails that might have spuriously high gradients
        % If we have an anomalously low or anomalously high value in profiles
        % they will result in high gradient but these might not be indicative
        % of background state
        Tsort(Tsort < prctile(Tvec, 5)) = nan;
        Tsort(Tsort > prctile(Tvec, 95)) = nan;

        mask = ~isnan(Tsort);
        Tmasked = Tsort(mask);
        zmasked = zthorpe(mask);
        [~,uinds] = unique(Tmasked);
        % uniform temperature!
        if length(uinds) == 1, continue; end

        % find location of chosen isotherms (Tbins) in sorted profiles
        zbins = interp1(Tmasked(uinds), zthorpe(uinds), Tbins);

        if optional
            zTj(:, ll) = interp1(Tmasked, zmasked, Tj);
        end

        % NaNs in zbins should only occur at the ends. If the profile didn't
        % contain some isotherms in Tbins then corresponding zbin is NaN...
        % We cannot use this profile to figure out distance between the
        % corresponding isopycnals
        dzmat(:, ll) = diff(zbins);

        if plotflag
            ct0 = find_approx(chi.time(chit0:chit1), vdisp.time(l0)) + chit0-1;
            ct1 = find_approx(chi.time(chit0:chit1), vdisp.time(l1)) + chit0-1;

            plot(htemp, Tvec, zvec, '-', 'linewidth', 0.5, 'HandleVisibility', 'off');
            hl = plot(hdisp, vdisp.time(l0:l1), zvec, '-', 'linewidth', 2);
            plot(hsort, Tsort, zthorpe, '-', 'HandleVisibility', 'off');
            plot(htemp2, vdisp.time(l0:l1), Tvec, '-', 'linewidth', 0.5);
            if ~chi_is_empty
                semilogy(hchi, chi.time(ct0:ct1), chi.chi(ct0:ct1), '-', ...
                         'handlevisibility', 'off');
            end
            if isnan(yloc)
                yloc = 0.95*min(vdisp.dis_z(t0:t1));
            end
            hsc = scatter(hchi_scatter, ...
                          chi.T(ct0:ct1), chi.chi(ct0:ct1), ...
                          32, repmat(hl.Color, [length(ct0:ct1), 1]), ...
                          'filled', 'HandleVisibility', 'off');

            tp_plot = Tp.tp(l0:l1);
            tp_plot(isnan(T.Tenh(l0:l1))) = nan;
            semilogy(htp, Tp.time(l0:l1), tp_plot, '-', 'color', hl.Color);
        end
    end

    % Now we reconstruct a mean "sorted profile" from binned dz/dT
    dz = nanmean(dzmat, 2)';
    numgood = sum(~isnan(dzmat), 2);
    % make sure the isotherm is present in at least three profiles
    dz(numgood < 3) = NaN;
    dzdT = dz./dT;

    wda.Tbins(1:length(Tbins), 1) = Tbins;
    wda.dTdz_bins(1:length(Tbins)-1, 1) = 1./dzdT;
    wda.dz(1:length(Tbins)-1, 1) = dz;

    if plotflag
        if ~chi_is_empty
            wda.dt = dt;
            % test that we can recover the value
            wda_proc = process_wda_estimate(chi, wda);
            % chiavg = isoscalar_average(chi.chi(chit0:chit1), chi.T(chit0:chit1), Tbins);
            % Jqavg = -chiavg .* dzdT * 4200 * 1025 * 0.5;
            % Ktavg = chiavg .* dzdT.^2 * 0.5;
            % Jqda = nansum(dz .* Jqavg)./nansum(dz);
            % if ~isnan(Jqda), assert(abs(Jqda - wda.Jq) < 1e-2), end
            title(hdisp, ['K_T = ' num2str(wda_proc.Kt, '%.1e') ...
                          ' | Jq = ' num2str(wda_proc.Jq)])
        end


        % plot all available gradients
        % fit T against z to get dT/dz == internal gradient
        % plot this as one reference
        [poly, ~, mu] = polyfit(zfull(~isnan(Tfull)), Tfull(~isnan(Tfull)), 1);
        Tzi = poly(1)/mu(2); % undo MATLAB scaling

        % our estimated mean "sorted" profile.
        dz(isnan(dz)) = 0;
        zprof = [mean(zthorpe)-std(zthorpe), mean(zthorpe)-std(zthorpe) + cumsum(dz)];

        for aa = [hsort, htemp]
            % internal
            plot(aa, polyval(poly, hsort.YLim, [], mu), aa.YLim, ...
                 'k-.', 'LineWidth', 2, 'DisplayName', ...
                 ['T_z^{fit} = ' num2str(Tzi, '%.1e') ' C/m']);

            % sorted
            plot(aa, Tbins, zprof, 'k-', 'linewidth', 2, 'displayname', ...
                 ['Mean T_z^{sort} = ' num2str(nanmean(1./dzdT), '%.1e') ' C/m']);

            % chi.dTdz
            plot(aa, ...
                 Tbins(1) + [0, diff(hsort.YLim) * mean(chi.dTdz(chit0:chit1))], ...
                 aa.YLim, 'k--', 'linewidth' ,2, 'DisplayName', ...
                 ['mean(chi.dTdz) = ' num2str(mean(chi.dTdz(chit0:chit1)), '%.1e')]);
        end

        % fit z against T to get dz/dT
        % [poly, ~, mu] = polyfit(Tfull, zfull, 1);
        % zTi = poly(1);
        % hline2 = plot(hsort, hsort.XLim, polyval(poly, hsort.XLim, [], mu), ...
        %               'k-', 'DisplayName', ['1/(dz/dT) = ' num2str(1./zTi, '%.1e')]);

        if optional
            % isoscalar average zTj : average over sorted profiles
            % this is usually junk.
            % I think the problem is choosing a sensible reference so that
            % we get a proper average.
            % instead I choose to average separations between isotherms
            % and then reconstruct a profile. That works a lot better.
            hline3 = plot(hsort, Tj, nanmean(zTj, 2), ...
                          'r-', 'linewidth', 2, 'DisplayName', ['W&DA']);
        end

        for aa=[hdisp, htemp2, htp]
            aa.XTickLabels = [];
            aa.XLabel.String = '';
        end

        hchi_scatter.XLabel.FontSize = fs * 1.1;

        % hsort.PlotBoxAspectRatio = [1 2 1];
        % htemp.PlotBoxAspectRatio = [1 2 1];
        hsort.XLim = [min(T.Tenh(t0:t1)) max(T.Tenh(t0:t1))];

        % hsort.Title.String = ['Jq_{DA} = ' num2str(Jqda, '%.1f') ...
        %                     ' | Jq_{i} = ' num2str(Jqi, '%.1f') ...
        %                     ' | Jq_{m} = ' num2str(Jqm, '%.1f')];
    end

end
function [data] = jitter(data, magnitude)
    data = data + magnitude * rand(size(data));
end

function [lims] = robust_lim(data)
    if max(abs(data)) > 3 * std(data)
        lims = [prctile(data, 1), prctile(data, 99)];
    else
        lims = [0.9 * min(data), 1.1 * max(data)];
    end
end