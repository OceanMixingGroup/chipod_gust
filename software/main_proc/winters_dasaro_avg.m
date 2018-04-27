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

    MIN_DIS_Z = 0.05; % minimum length (in metres) of a single up- or down-cast
    nquantiles = round(5/60 * dt); % effectively number of bins

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

    % Lets make some T bins by splitting into quantiles so that each bin has equal number of observations
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

    if plotflag
        CreateFigure;
        htemp = subplot(3,4,[1, 5, 9]); hold on; xlabel('Temp'); ylabel('z (accel)')
        hsort = subplot(3,4,[4, 8, 12]); hold on; xlabel('Temp (sorted)'); ylabel('z (interp)')
        hsort.XGrid = 'on'; hsort.XTick = Tbins;
        hsort.XAxis.TickLabelRotation = 30; hsort.GridAlpha = 1;
        hdisp = subplot(3,4,[2, 3]); hold on; xlabel('time'); ylabel('Displacement')
        hdisp.Color = 'none';
        if ~chi_is_empty
            hchi = subplot(3,4,[6,7]); hold on; xlabel('time'); ylabel(['\chi (1 sec, color)'])
        end
        htemp2 = axes('Position', hdisp.Position, 'Color', 'none'); hold on;
        htp = subplot(3,4,[10, 11]); hold on; xlabel(['time ' datestr(chi.time(chit0))]); ylabel(['Tp'])

        yloc = NaN;
    end

    for ll=1:length(lstart) % loop over "profiles"
        l0 = lstart(ll); l1 = lstop(ll);

        if abs(vdisp.dis_z(l0) - vdisp.dis_z(l1)) < MIN_DIS_Z, continue; end

        zvec = -vdisp.dis_z(l0:l1);
        Tvec = T.Tenh(l0:l1); %T.T(l0:l1);

        if all(isnan(Tvec)), continue; end

        if any(Tp.tp(l0:l1) > 15), continue; end

        % interpolate to uniform depths before sorting
        Tinterp = interp1(zvec, Tvec, zthorpe);

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

            plot(htemp, Tvec, zvec, '-', 'linewidth', 0.5);
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
            hsc = scatter(hsort, chi.T(ct0:ct1), yloc * ones(size(chi.T(ct0:ct1))), ...
                          200, hl.Color, 'HandleVisibility', 'off');

            tp_plot = Tp.tp(l0:l1);
            tp_plot(isnan(T.Tenh(l0:l1))) = nan;
            semilogy(htp, Tp.time(l0:l1), tp_plot, '-', 'color', hl.Color);
        end
    end

    % reconstruct profile from binned dz/dT
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

        % fit T against z to get dT/dz == internal gradient
        [poly, ~, mu] = polyfit(zfull(~isnan(Tfull)), Tfull(~isnan(Tfull)), 1);
        Tzi = poly(1)/mu(2); % undo MATLAB scaling
        if plotflag
            hline = plot(hsort, polyval(poly, hsort.YLim, [], mu), hsort.YLim, ...
                         'k-', 'DisplayName', ['dTdz_i = ' num2str(Tzi, '%.1e')]);
            plot(htemp, polyval(poly, hsort.YLim, [], mu), htemp.YLim, ...
                 'k-', 'DisplayName', ['dTdz_i = ' num2str(Tzi, '%.1e')]);
        end

        % fit z against T to get dz/dT
        % [poly, ~, mu] = polyfit(Tfull, zfull, 1);
        % zTi = poly(1);
        % hline2 = plot(hsort, hsort.XLim, polyval(poly, hsort.XLim, [], mu), ...
        %               'k-', 'DisplayName', ['1/(dz/dT) = ' num2str(1./zTi, '%.1e')]);

        if optional
            % isoscalar average zTj
            hline3 = plot(hsort, Tj, nanmean(zTj, 2), ...
                          'r-', 'linewidth', 2, 'DisplayName', ['W&DA']);
        end

        dz(isnan(dz)) = 0;
        zprof = [mean(zthorpe)-std(zthorpe), mean(zthorpe)-std(zthorpe) + cumsum(dz)];
        hline4 = plot(hsort, Tbins, zprof, 'k-', 'linewidth', 2, ...
                      'displayname', ['average \Delta z mean=' ...
                            num2str(nanmean(1./dzdT), '%.1e')]);

        if isfield(chi, 'dTdz')
            hline5 = plot(hsort, Tbins(1) + [0, diff(hsort.YLim) * mean(chi.dTdz(chit0:chit1))], ...
                          hsort.YLim, 'r-', 'linewidth' ,2, 'DisplayName', ...
                          ['mean(chi.dTdz) = ' num2str(mean(chi.dTdz(chit0:chit1)), '%.1e')]);
        end

        hzfull = plot(hdisp, vdisp.time(t0:t1), ...
                      -(vdisp.dis_z(t0:t1) - vdisp.dis_z(t0)), ...
                      'color', [1 1 1]*0.75, 'linewidth', 2);
        if ~chi_is_empty
            hcfull = plot(hchi, chi.time(chit0:chit1), ...
                          chi.chi(chit0:chit1), 'HandleVisibility', 'off', ...
                          'color', [1 1 1]*0.75, 'linewidth', 2);
            uistack(hcfull, 'bottom')
        end
        uistack(hzfull, 'bottom')

        hscfull = scatter(hsort, Tchi, yloc * ones(size(chi.T(chit0:chit1))), ...
                          200, 'k', 'HandleVisibility', 'off');

        hfullT = plot(htemp2, T.time, T.T, 'color', [1 1 1]*0.5);
        ylim(htemp2, robust_lim(chi.T(chit0:chit1)))
        uistack(hfullT, 'bottom');

        hchi.YScale = 'log';
        hcloud = plot(hsort, Tfull, zfull, '.', ...
                      'color', [1 1 1]*0.8, 'HandleVisibility', 'off');
        uistack(hcloud, 'bottom');

        tpt0 = find_approx(Tp.time, vdisp.time(t0));
        tpt1 = find_approx(Tp.time, vdisp.time(t1));
        htpfull = semilogy(htp, Tp.time(tpt0:tpt1), Tp.tp(tpt0:tpt1), 'color', ...
                           [1 1 1]*0.6);
        ylim(htp, robust_lim(Tp.tp(tpt0:tpt1)))
        uistack(htpfull, 'bottom')
        axes(htp); datetick('x', 'MM:SS', 'keeplimits');

        if ~chi_is_empty
            semilogy(hchi, chi.time(chit0:chit1), chi.spec_area(chit0:chit1), ...
                     'color', 'k', 'linewidth', 2, 'displayname', ['spec_area'])
            legend(hchi, '-dynamiclegend')
            plot(hchi, hchi.XLim, [1, 1]*chi.spec_floor * chi.nfft, 'k-', ...
                 'displayname', 'noise floor');
            plot(hchi, hchi.XLim, 2*[1, 1]*chi.spec_floor * chi.nfft, 'k--', ...
                 'displayname', '2x noise floor');
            axes(hchi); datetick('x', 'MM:SS', 'keeplimits');
        end

        hleg = legend(hsort, '-dynamiclegend');
        hleg.Position = [0.7479    0.7733    0.1650    0.0800];
        axes(hdisp); datetick('x', 'MM:SS', 'keeplimits');
        uistack(htemp2, 'bottom');
        if ~chi_is_empty
            linkaxes([hchi hdisp htp htemp2], 'x');
        else
            linkaxes([hdisp htp htemp2], 'x');
        end

        hdisp.XLim = chi.time([chit0, chit1]);
        htemp2.YAxisLocation = 'right';
        htemp2.XTickLabels = [];
        htemp2.XTick = hdisp.XTick;

        hsort.PlotBoxAspectRatio = [1 2 1];
        htemp.PlotBoxAspectRatio = [1 2 1];
        hsort.XLim = [min(T.Tenh(t0:t1))-1e-3 max(T.Tenh(t0:t1))+1e-3];

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
        lims = [prctile(data, 2), prctile(data, 98)];
    else
        lims = [0.9 * min(data), 1.1 * max(data)];
    end
end