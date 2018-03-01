% Given appropriate temperature gradients in wda, estimate volume avg Kt, Jq
function [wda] = process_wda_estimate(chi, wda)
    ticwda = tic;
    disp('Processing Winters & D''Asaro estimate. Takes 120s for 1 year.');

    ndt = round(wda.dt(1) / diff(chi.time(1:2)*86400));

    % at this point, chi might have been trimmed (it is trimmed in
    % combine_turbulence); so we need to find the right time chunks in wda so
    % that things line up right.
    % chunk0 is the wda index of the time chunk that is closest to the first good chi observation
    chunk0 = find(wda.tstart > chi.time(1), 1, 'first');
    chit0 = find_approx(chi.time, wda.tstart(chunk0));

    wda.Jq = nan(1, length(wda.tstart));
    wda.Kt = nan(1, length(wda.tstart));
    wda.eps = nan(1, length(wda.tstart));
    wda.chi = nan(1, length(wda.tstart));
    wda.dTdz = nan(1, length(wda.tstart));

    for tt=chunk0:length(wda.tstart)
        maxt = min(chit0+2*ndt, length(chi.time));

        % narrow the search space for speed using ndt
        chit0 = find_approx(chi.time(chit0:maxt), wda.tstart(tt)) + chit0-1;
        chit1 = find_approx(chi.time(chit0+1:maxt), wda.tstop(tt), 2) + chit0;

        % if at the end, they probably don't line up; discard this one point
        if chit1(1) == length(chi.chi), break; end

        if abs(chi.time(chit1(1)) - wda.tstop(tt))*86400 > 0.5
            chit1 = chit1(2);
        else
            chit1 = chit1(1);
        end

        % make sure the time chunks are lining up
        % disp(['Compare chit0 = ' datestr(chi.time(chit0)) ' vs tstart = ' datestr(wda.tstart(tt))]);
        % disp(['Compare chit1 = ' datestr(chi.time(chit1)) ' vs tstop = ' datestr(wda.tstop(tt))]);
        assert(abs(chi.time(chit0) - wda.tstart(tt))*86400 < 0.9)
        assert(abs(chi.time(chit1) - wda.tstop(tt))*86400 < 0.9)

        if all(isnan(chi.chi(chit0:chit1))), continue; end

        Tchi = chi.T(chit0:chit1);

        Tbins = wda.Tbins(:, tt);
        Tbins = Tbins(~isnan(Tbins));

        if length(Tbins) > 1
            chiavg = isoscalar_average(chi.chi(chit0:chit1), Tchi, Tbins);
            epsavg = isoscalar_average(chi.eps(chit0:chit1), Tchi, Tbins);

            dz = wda.dz(1:length(Tbins)-1, tt);
            dTdz = wda.dTdz_bins(1:length(Tbins)-1, tt);
            Jqvec = -chiavg' ./ dTdz * 4200 * 1025 * 0.5;
            Ktvec = chiavg' ./ dTdz.^2 * 0.5;

            wda.chi(tt) = sum(dz .* chiavg', 'omitnan')./sum(dz, 'omitnan');
            wda.eps(tt) = sum(dz .* epsavg', 'omitnan')./sum(dz, 'omitnan');
            wda.Jq(tt) = sum(dz .* Jqvec, 'omitnan')./sum(dz, 'omitnan');
            wda.Kt(tt) = sum(dz .* Ktvec, 'omitnan')./sum(dz, 'omitnan');
            wda.dTdz(tt) = sum(dz .* dTdz, 'omitnan')./sum(dz, 'omitnan');
        end
    end

    wda.time = (wda.tstart + wda.tstop)/2;
    toc(ticwda);
end