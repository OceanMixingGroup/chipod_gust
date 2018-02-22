% [wda] = winters_dasaro_avg(t0, t1, data, chi, T, plotflag)
% Treats the chipod as a profiler and applies the method of Winters & D'Asaro to estimate heat flux.
% Inputs:
%        t0, t1 - time indices for accelerometer data (usually 50Hz)
%        data, chi, T - structures passed in from chi_main_driver
%        plotflag - if 1, make informative plot showing resorted profile.
% Outputs:
%        wda - structure with fields
%           tstart, tstop - start,end of time chunk
%

function [wda] = winters_dasaro_avg(t0, t1, data, chi, T, plotflag)

    optional = 0;

    tstart = data.time(t0);
    tstop = data.time(t1);

    zfull = -(data.a_dis_z(t0:t1)-data.a_dis_z(t0));
    Tfull = T.T(t0:t1);

    % determine "profile" boundaries
    [~, locs1] = findpeaks(data.a_dis_z(t0:t1));
    [~, locs2] = findpeaks(-data.a_dis_z(t0:t1));
    locs = sort([t0; locs1+t0-1; locs2+t0-1; t1]);

    chit0 = find_approx(chi.time, tstart);
    chit1 = find_approx(chi.time, tstop);
    Tchi = chi.T(chit0:chit1);

    il = 1; Tprof = [];
    for ll=1:length(locs)-1
        l0 = locs(ll); l1 = locs(ll+1);
        if abs(data.a_dis_z(l0) - data.a_dis_z(l1)) > 0.1
            lstart(il) = locs(ll);
            lstop(il) = locs(ll+1);
            il = il+1;

            ct0 = find_approx(chi.time(chit0:chit1), data.time(l0)) + chit0-1;
            ct1 = find_approx(chi.time(chit0:chit1), data.time(l1)) + chit0-1;

            % save 1sec avg temp in valid profiles. This is used to make bins
            Tprof = [Tprof chi.T(ct0:ct1)];
        end
    end

    % Lets make some T bins by splitting into quantiles so that each bin has equal number of observations
    nquantiles = 5;
    Tbins = generate_wda_bins(Tprof, nquantiles);
    dT = diff(Tbins);

    % initialized returned structure
    wda.nbins = nquantiles+2;
    wda.Tbins = (nan(nquantiles+2, 1));
    wda.tstart = tstart;
    wda.tstop = tstop;
    wda.dTdz = (nan(nquantiles+2, 1));
    wda.dz = (nan(nquantiles+2, 1));

    % if no long enough profiles are found
    if isempty(Tprof), return; end

    zthorpe = linspace(min(zfull), max(zfull), 1000);
    if optional
        Tj = linspace(min(Tfull), max(Tfull), 1000);
        zTj = nan(length(Tj), length(locs));
    end
    dzmat = nan(length(dT), length(locs));

    if plotflag
        clf;
        htemp = subplot(2,4,[1, 5]); hold on; xlabel('Temp'); ylabel('z (accel)')
        hsort = subplot(2,4,[4, 8]); hold on; xlabel('Temp (sorted)'); ylabel('z (interp)')
        hsort.XGrid = 'on'; hsort.XTick = Tbins; hsort.XAxis.TickLabelRotation = 30;
        hdisp = subplot(2,4,[2, 3]); hold on; xlabel('time'); ylabel('Displacement')
        hchi = subplot(2,4,[6,7]); hold on; xlabel('time'); ylabel(['\chi (1 sec)'])
        htemp2 = axes('Position', hdisp.Position, 'Color', 'none'); hold on;
        hdisp.Color = 'none';
    end

    yloc = NaN;
    for ll=1:length(lstart) % loop over "profiles"
        l0 = lstart(ll); l1 = lstop(ll);

        zvec = -(data.a_dis_z(l0:l1) - data.a_dis_z(locs(1)));
        Tvec = T.T(l0:l1);

        % interpolate to uniform depths before sorting
        % Tsort is on depth grid zthorpe
        Tinterp = interp1(zvec, Tvec, zthorpe);
        [Tsort, inds] = thorpeSort(Tinterp);

        mask = ~isnan(inds);
        % add some jitter so the interpolation works when gradients are low
        % save these for averaging over profiles later
        if optional
            zTj(:, ll) = interp1(jitter(Tsort(mask), 1e-5), zthorpe(mask), Tj);
        end

        try
            zbins = interp1(jitter(Tsort(mask), 1e-5), zthorpe(mask), Tbins);
        catch ME
            % if jitter didn't work by chance try again
            zbins = interp1(jitter(Tsort(mask), 1e-5), zthorpe(mask), Tbins);
        end
        % NaNs in zbins should only occur at the ends. If the profile didn't
        % contain some isotherms in Tbins then corresponding zbin is NaN...
        % We cannot use this profile to figure out distance between the
        % corresponding isopycnals
        dzmat(:, ll) = diff(zbins);

        if plotflag
            ct0 = find_approx(chi.time(chit0:chit1), data.time(l0)) + chit0-1;
            ct1 = find_approx(chi.time(chit0:chit1), data.time(l1)) + chit0-1;

            plot(htemp, Tvec, zvec, '-', 'linewidth', 0.5);
            hl = plot(hdisp, data.time(l0:l1), zvec, '-', 'linewidth', 2);
            plot(hsort, Tsort, zthorpe, '-', 'HandleVisibility', 'off');
            plot(htemp2, data.time(l0:l1), Tvec, '-', 'linewidth', 0.5);
            semilogy(hchi, chi.time(ct0:ct1), chi.chi(ct0:ct1), '-');
            if isnan(yloc)
                yloc = 0.95*min(data.a_dis_z(t0:t1));
            end
            hsc = scatter(hsort, chi.T(ct0:ct1), yloc * ones(size(chi.T(ct0:ct1))), ...
                          200, hl.Color, 'HandleVisibility', 'off');

        end
    end

    % reconstruct profile from binned dz/dT
    dz = nanmean(dzmat, 2)';
    numgood = sum(~isnan(dzmat), 2);
    dz(numgood < 3) = NaN;
    dzdT = dz./dT;

    wda.Tbins(1:length(Tbins), 1) = Tbins;
    wda.dTdz(1:length(Tbins)-1, 1) = 1./dzdT;
    wda.dz(1:length(Tbins)-1, 1) = dz;

    % test that we can recover the value
    % wda = process_wda_estimate(chi, wda);
    % chiavg = isoscalar_average(chi.chi(chit0:chit1), chi.T(chit0:chit1), Tbins);
    % Jqavg = -chiavg .* dzdT * 4200 * 1025 * 0.5;
    % Ktavg = chiavg .* dzdT.^2 * 0.5;
    % Jqda = nansum(dz .* Jqavg)./nansum(dz);
    % if ~isnan(Jqda), assert(abs(Jqda - wda.Jq) < 1e-2), end

    if plotflag
        % fit T against z to get dT/dz
        [poly, ~, mu] = polyfit(zfull, Tfull, 1);
        Tzi = poly(1)/mu(2); % undo MATLAB scaling
        if plotflag
            hline = plot(hsort, polyval(poly, hsort.YLim, [], mu), hsort.YLim, ...
                         'k-', 'DisplayName', ['dTdz_i = ' num2str(Tzi, '%.1e')]);
            plot(htemp, polyval(poly, hsort.YLim, [], mu), htemp.YLim, ...
                 'k-', 'DisplayName', ['dTdz_i = ' num2str(Tzi, '%.1e')]);
        end

        % fit z against T to get dz/dT
        [poly, ~, mu] = polyfit(Tfull, zfull, 1);
        zTi = poly(1);
        hline2 = plot(hsort, hsort.XLim, polyval(poly, hsort.XLim, [], mu), ...
                      'k-', 'DisplayName', ['1/(dz/dT) = ' num2str(1./zTi, '%.1e')]);

        if optional
            % isoscalar average zTj
            hline3 = plot(hsort, Tj, nanmean(zTj, 2), ...
                          'r-', 'linewidth', 2, 'DisplayName', ['W&DA']);
        end

        dz(isnan(dz)) = 0;
        zprof = [mean(zthorpe)-std(zthorpe), mean(zthorpe)-std(zthorpe) + cumsum(dz)];
        hline4 = plot(hsort, Tbins, zprof, 'k-', 'linewidth', 2, ...
                      'displayname', ['average \Delta z mdn=' ...
                            num2str(nanmean(1./dzdT), '%.1e')]);

        hzfull = plot(hdisp, data.time(t0:t1), ...
                      -(data.a_dis_z(t0:t1) - data.a_dis_z(t0)), ...
                      'color', [1 1 1]*0.75, 'linewidth', 2);
        hcfull = plot(hchi, chi.time(chit0:chit1), ...
                      chi.chi(chit0:chit1), ...
                      'color', [1 1 1]*0.75, 'linewidth', 2);
        uistack(hzfull, 'bottom')
        uistack(hcfull, 'bottom')

        hscfull = scatter(hsort, Tchi, yloc * ones(size(chi.T(chit0:chit1))), ...
                          200, 'k', 'HandleVisibility', 'off');

        hchi.YScale = 'log';
        hcloud = plot(hsort, Tfull, zfull, '.', ...
                      'color', [1 1 1]*0.8, 'HandleVisibility', 'off');
        uistack(hcloud, 'bottom');

        legend(hsort, '-dynamiclegend', 'Location', 'NorthWest');
        axes(hchi); datetick('x', 'MM:SS', 'keeplimits');
        axes(hdisp); datetick('x', 'MM:SS', 'keeplimits');
        uistack(htemp2, 'bottom');
        linkaxes([hchi hdisp htemp2], 'x');
        hdisp.XLim = chi.time([chit0, chit1]);
        htemp2.YAxisLocation = 'right';
        htemp2.XTickLabels = [];
        htemp2.XTick = hdisp.XTick;

        hsort.PlotBoxAspectRatio = [1 2 1];
        htemp.PlotBoxAspectRatio = [1 2 1];
        hsort.XLim = [min(T.T(t0:t1))-1e-3 max(T.T(t0:t1))+1e-3];
        hsort.Title.String = ['Jq_{DA} = ' num2str(Jqda, '%.1f') ...
                            ' | Jq_{i} = ' num2str(Jqi, '%.1f') ...
                            ' | Jq_{m} = ' num2str(Jqm, '%.1f')];
    end
end

function [data] = jitter(data, magnitude)
    data = data + magnitude * rand(size(data));
end