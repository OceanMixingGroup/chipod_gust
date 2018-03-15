% [Tp_floor, T_floor] = get_nosie_floors(Tvolt, Tpcoef)
%  Input:
%         Tvolt: Raw temperature voltage time series
%         Tcoef: Temperature calibration coefficient
%         Tpcoef: Tp calibration coefficients
%       Fs: sampling frequency
%
%  Output:
%         Tp_floor: spectral level floor
%         T_floor : temperature signal floor (C)

function [Tp_floor, T_floor] = get_noise_floors(Tvolt, Tcoef, Tpcoef, Fs, ...
                                                       Tvec, Tpvec)

    if Tpcoef(2) < 1e-5
        gain = Tpcoef(1); % gusT
    else
        gain = Tpcoef(2); % chipod
    end

    assert(gain < 0.5, 'get_tp_spec_floor: Tp gain seems wrong.')

    % Becherer & Moum (2017) : eqns 27-31
    % We need to scale their CTp by TP gain i.e. Tpcoef(1) or Tpcoef(2)
    % (see data_reduction stuff)

    ctp = (2 * nanmean(Tvolt) * Tcoef(3) + Tcoef(2))/gain;

    % 4.096 V / (16 bit quantization)
    % + A2D error of 0-6 counts peak-peak == 3 counts amplitude
    Tp_floor = ctp.^2 * ((4.096/2^16) * 3)^2 / Fs;

    T_floor = ctp * 4.096/2^16 * gain;

    dT = abs(diff(Tvec));
    dTmin = nanmin(dT(dT > 0));
    if abs(T_floor - dTmin)/dTmin * 100 > 1
        disp(['Warning: estimated temperature noise floor is different from ' ...
              'min(abs(diff(T)) by more than 1% '])
    end

    % generate bit noise time series, and figure out variance when
    % differentiator hits noise floor
    % 4.096 V / (16 bit quantization) + A2D error of 0-6 counts peak-peak == 3 counts amplitude
    bit_noise_volts = randi([0, 6], length(Tvolt), 1) * (4.096/2^16);

    tpnoise = calibrate_tp(bit_noise_volts, Tpcoef, Tvolt, Tcoef, 100*ones(size(Tpvec)) );

    ndt = Fs;
    floor1 = moving_var(tpnoise, ndt, ndt);
    est_TP_floor = prctile(floor1, 95);

    % integral (noise_floor * df) = variance
    % noise_floor * integral_(0)^(50 Hz) (df) = variance
    if abs(Tp_floor * Fs/2 - est_TP_floor)/est_TP_floor * 100 > 5
        disp(['Warning: estimated Tp noise floor is different from ' ...
              'bit noise variance by more than 5% '])
    end

    % noise floor debugging plots
    % [S, f] = fast_psd(tpnoise(1:10*Fs), Fs/2, Fs)
    % loglog(f, S); hold on;
    % % these two should be roughly mean of estimated spectrum
    % plot(xlim, [1, 1] * est_TP_floor/Fs*2, 'k-')
    % plot(xlim, [1, 1] * Tp_floor, 'k-')

    % mask = moving_var(chi.T2Pt, ndt, 1) < 1.2*(chi.T2P_spec_floor * 50);
    % t2v = moving_var(chi.T2Pt, ndt, ndt);
    % figure;
    % ax(1) = subplot(311);
    % semilogy(chi.time_tp(1:ndt:end), t2v); hold on;
    % semilogy(chi.time_tp(1:ndt:end), floor2);
    % liney(chi.T2P_spec_floor * 50)
    % ax(2) = subplot(312); plot(chi.time, chi.T2);
    % ax(3) = subplot(313); plot(chi.time_tp, chi.T2Pt); hold on; ...
    %         plot(chi.time_tp(mask), chi.T2Pt(mask), '.')
    % linkaxes(ax, 'x')
end