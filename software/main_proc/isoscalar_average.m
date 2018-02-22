function [avg, bins] = isoscalar_average(data, scalar, bins)

    if isscalar(bins)
        bins = linspace(min(scalar), max(scalar), bins);
    elseif isvector(bins)
        bins = bins;
    end

    Y = discretize(scalar, bins, 'IncludedEdge', 'left');

    avg = nan(1, length(bins)-1);
    for binindex=1:length(bins)-1
        if sum(~isnan(data(Y == binindex))) > 2
            avg(binindex) = nanmean(data(Y == binindex));
        end
    end
end