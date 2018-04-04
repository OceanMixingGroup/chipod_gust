% Estimate a sign for the sorted dTdz as the sign(hourly median(reference dTdz))
%
%         function [sgn] = get_wda_sign(time_vec, refTz, min_dTdz)
%
% Inputs:
%       time_vec: time vector for the sorted gradient
%       Tz : reference dTdz to extract sign from
%       min_dTdz : when less than minimum value, use a 2 hour median
% Outputs:
%       sgn : vector of +/-1 indicating sign

function [sgn] = get_wda_sign(time_vec, refTz, min_dTdz)

    if ~exist('min_dTdz', 'var'), min_dTdz = 1e-3; end

    dt = round((refTz.time(2)-refTz.time(1))*86400);
    refTzi = interp1(refTz.time, refTz.Tz, time_vec);

    % sign of hourly moving median
    sgn = sign(movmedian(refTzi, 60*60/dt, 'omitnan'));
    sgn2h = sign(movmedian(refTzi, 2*60*60/dt, 'omitnan'));

    % If T_z is crossing CP.min_dTdz, we shouldn't take it too seriously
    Tz_min_cross = generate_min_dTdz_crossing_mask(refTzi, min_dTdz, 0);
    sgn(Tz_min_cross) = nan;
    sgn = fillmissing(sgn, 'nearest'); % nearest-neighbour filling only works for tiny gaps
    sgn(abs(refTzi) < min_dTdz) = nan;
    sgn(isnan(sgn)) = sgn2h(isnan(sgn));
end
