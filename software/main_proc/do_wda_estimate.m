function [chi_wda] = do_wda_estimate(params, data, chi, T, Tp)

    % ENHANCE! T and Tp to get higher resolution time series
    % Mudge & Lueck (1994)
    % cross-over frequency at 10Hz

    nanT = isnan(T.T);
    nanTp = isnan(Tp.tp);

    Tenh = nan(size(Tp.tp));

    if all(nanT) || all(nanTp)
        disp(['do_wda_estimate received all-NaN T or all-NaN Tp']);
    else
        T.T(~nanT) = fillmissing(deglitch(T.T(~nanT), 100, 5, 'b', 1), 'linear');
        Tinterp = fillmissing(interp1(T.time, T.T, Tp.time, 'nearest'), 'nearest');
        Tenh(~nanTp) = combine_ttc(Tinterp(~nanTp), Tp.tp(~nanTp), 2*pi*10, 100, 0)';
    end

    % % find out when T sensor is at noise floor
    % mask = movvar(T.T, 50*2.5) < T.floor^2 * 9/7;
    % Tmask = interp1(T.time, double(mask), Tp.time, 'nearest');
    % Tmask(isnan(Tmask)) = 0;
    % Tmask = logical(Tmask);

    % % find where there is no detectable temperature fluctuation
    % mask = chi.spec_area < 2*chi.spec_floor*chi.nfft;
    % tpmask = interp1(chi.time, double(mask), Tp.time, 'nearest');
    % tpmask(isnan(tpmask)) = 0;
    % tpmask = logical(tpmask);

    Temasked = Tenh;
    % plot(Tp.time, Tenh); hold on;
    % when Tp is at noise floor, replace Tenh with 50Hz T time series (interpolated)
    % Temasked(tpmask) = Tinterp(tpmask);
    % % plot(Tp.time, Temasked);
    % % when both T, Tp are at noise floors; T is constant as far as we can tell.
    % Temasked(Tmask & tpmask) = nan;
    % Temasked = fillmissing(Temasked, 'previous');
    % % plot(Tp.time, Temasked, 'k')

    % use enhanced temperature to get chi.T
    chi.T = moving_average(Temasked', 100, 100);
    T.Tenh = Temasked;
    T.time_enh = Tp.time;

    % need accelerometer/pressure data at 100Hz too since T has been enhanced
    vdisp.time = T.time_enh;
    if params.do_P
        vdisp.dis_z = interp1(data.time, data.p_dis_z, T.time_enh);
    else
        vdisp.dis_z = interp1(data.time, data.a_dis_z, T.time_enh);
    end

    ndt = params.wda_dt * round(1/diff(vdisp.time(1:2)*86400));
    idx = 1;
    plotflag = 0;
    Nt = length(1:ndt:length(vdisp.dis_z));
    wda = cell(Nt, 1);
    for t0=1:ndt:length(vdisp.dis_z)
        wda{idx} = winters_dasaro_avg(t0, min(t0 + ndt, length(vdisp.dis_z)), ...
                                      vdisp, chi, T, Tp, params.wda_dt, ...
                                      plotflag);
        idx = idx+1;
    end

    chi_wda = merge_cell_structs(wda, 1);
    chi_wda.dt = params.wda_dt;
    chi_wda.nbins = chi_wda.nbins(1);

    % keyboard;
    % wda_proc = process_wda_estimate(chi, chi_wda);
    % plot_estimate(wda_proc)

    % old 50Hz style inference
    % acc50hz.time = data.time;
    % acc50hz.dis_z = data.a_dis_z;
    % T50hz.Tenh = T.T;
    % T50hz.time_enh = T.time;
    % wda50 = do_wda_estimate(params.wda_dt, acc50hz, chi, T50hz, Tp);

    % wda_proc100 = process_wda_estimate(chi, chi.wda);
    % wda_proc50 = process_wda_estimate(chi, wda50);

    % load ../proc/T_m.mat

    % CreateFigure;
    % ax1 = subplot(311);
    % plot(T.time, T.T); hold on;
    % plot(Tp.time, Tenh, 'k');
    % plot(T1.time, T1.T)
    % plot(T2.time, T2.T)
    % legend('50Hz T', 'Enhanced 100Hz T', 'CTD1', 'CTD2')
    % datetick('x', 'mm/dd HH:MM', 'keeplimits')

    % ax2 = subplot(312); hold on;
    % plot(wda_proc50.time, wda_proc50.dTdz);
    % plot(wda_proc100.time, wda_proc100.dTdz);
    % % plot(wda_proc2.time, wda_proc2.dTdz);
    % % plot(wda_proc3.time, wda_proc3.dTdz, 'k');
    % plot(chi.time, abs(chi.dTdz), 'k')
    % title('different dT/dz estimates')
    % legend('old 50Hz', 'pre-emphasized 100Hz', '|moor|');
    % datetick('x', 'mm/dd HH:MM', 'keeplimits')

    % ax3 = subplot(313); hold on;
    % plot(Tp.time, Tp.tp);
    % plot(Tp.time, tpsmall, 'color', [1 1 1]*0.65);

    % % plot(wda_proc2.time, wda_proc2.dTdz);
    % % plot(wda_proc3.time, wda_proc3.dTdz, 'k');
    % % plot(chi.time, abs(chi.dTdz), 'k')
    % datetick('x', 'mm/dd HH:MM', 'keeplimits')
    % linkaxes([ax1, ax2, ax3], 'x');
    % xlim([Tp.time(1) Tp.time(end)])

    % keyboard;
    % savefig(['~/rama/images/526-' rfid(5:10) '.fig']);
    % keyboard;
    % wda_compare_plot; % script to make comparison plots
end