function [Tbins] = generate_wda_bins(Tchi, nquantiles)

    % temperature resolution
    % 1. If I use min/max as edges then because we only have one observation of max/min, that bin is useless.
    %    So I adjust based on dTres to make better use of data. (really minor tweak)
    % 2. If bin spacing is less than that; lets remove those bin edges
    % because we work file-by-file, can't really determine this empirically
    dTres = 5e-5; % 1mK precision

    % adjust min/max edges so that we use as much data as possible
    Tbins = (sort([prctile(Tchi, 5), quantile(Tchi, nquantiles), prctile(Tchi, 95)]));

    % make sure bin widths are > dTres
    bins2 = Tbins;
    ii = 1;
    while ii <= length(Tbins)-1
        for jj=ii+1:length(Tbins)
            if Tbins(jj) - Tbins(ii) < dTres
                bins2(jj) = nan;
            else
                break;
            end
        end
        ii = jj;
    end

    Tbins = bins2(~isnan(bins2));
end