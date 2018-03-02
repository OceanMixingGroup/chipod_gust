function [Tz_min_cross] = generate_min_dTdz_crossing_mask(dTdz, min_dTdz, debug_plot)
% [Tz_min_cross] = generate_min_dTdz_crossing_mask(dTdz, min_dTdz, debug_plot)
%
% Generates a mask for points that are on either side of a crossing of
% min_dTdz. I.e., if Tz(t0) < 1e-3 and Tz(t0-1), Tz(t0+1) > 1e-3, we might
% not want to keep t0-1, t0+1 values.
%
% Inputs:
%         dTdz - dTdz timeseries
%     min_dTdz - threshold
%   debug_plot - if 1, make debugging plot

    if size(dTdz, 1) > 1
        transpose = 1;
        dTdz = dTdz';
    else
        transpose = 0;
    end

    sgn = sign(abs(dTdz) - min_dTdz) .* sign(dTdz);
    sgnflip = flip(sgn, 2);
    sTz = [1 sgn(1:end-1) .* sgn(2:end)];
    sTz1 = flip([1 sgnflip(1:end-1) .* sgnflip(2:end)], 2);

    % debugging plots
    if debug_plot
        Tzmcross = ((sTz == -1) | (sTz1 == -1) | (abs(dTdz) < min_dTdz));
        Tz2 = dTdz; Tz2(Tzmcross) = nan;

        figure;
        plot(dTdz);
        hold on;
        plot(Tz2, 'b-')
        liney([-min_dTdz, min_dTdz, 0])
    end

    Tz_min_cross = ((sTz == -1) | (sTz1 == -1));

    if transpose, Tz_min_cross = Tz_min_cross'; end
end