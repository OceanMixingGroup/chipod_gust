% Given appropriate temperature gradients in wda, estimate volume avg Kt, Jq
function [wda] = process_wda_estimate(chi, wda)

    for tt=1:length(wda.tstart)
        chit0 = find_approx(chi.time, wda.tstart(tt));
        chit1 = find_approx(chi.time, wda.tstop(tt));
        Tchi = chi.T(chit0:chit1);

        Tbins = wda.Tbins(:, tt);
        Tbins = Tbins(~isnan(Tbins));

        chiavg = isoscalar_average(chi.chi(chit0:chit1), Tchi, Tbins);

        dz = wda.dz(1:length(Tbins)-1, tt);
        dTdz = wda.dTdz(1:length(Tbins)-1, tt);
        Jqvec = -chiavg' ./ dTdz * 4200 * 1025 * 0.5;
        Ktvec = chiavg' ./ dTdz.^2 * 0.5;

        wda.Jq(tt) = nansum(dz .* Jqvec)./nansum(dz);
        wda.Kt(tt) = nansum(dz .* Ktvec)./nansum(dz);
    end

    wda.time = (wda.tstart + wda.tstop)/2;
end