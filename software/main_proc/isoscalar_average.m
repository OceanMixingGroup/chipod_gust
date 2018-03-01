function [avg, bins] = isoscalar_average(data, scalar, bins)

    if isscalar(bins)
        bins = linspace(min(scalar), max(scalar), bins);
    elseif isvector(bins)
        bins = bins;
    end

    Y = discretize(scalar, bins);
    avg = nan(1, length(bins)-1);

    for binindex=1:length(bins)-1
        databin = data(Y == binindex);
        numvalid = sum(~isnan(databin));
        if numvalid > 2
            avg(binindex) = sum(databin, 'omitnan')/numvalid;
        end
    end
end