function [avg, bins] = isoscalar_average(data, scalar, bins)

    if isscalar(bins)
        bins = linspace(min(scalar), max(scalar), bins);
    elseif isvector(bins)
        bins = bins;
    end

    Y = discretize(scalar, bins); % this is the slow step
    avg = nan(size(data, 1), length(bins)-1);

    for binindex=1:length(bins)-1
        databin = data(:, Y == binindex);
        numvalid = sum(~isnan(databin), 2);
        if numvalid(1) > 2
            avg(:, binindex) = sum(databin, 2, 'omitnan')/numvalid(1);
        end
    end
end