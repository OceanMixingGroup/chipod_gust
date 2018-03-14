function [chi, percentage] = ApplyMask(chi, maskvar, relation, criterion, ...
                                       maskname, Trange, fill_value)

    if ~exist('Trange', 'var') | isempty(Trange), Trange = 1:length(chi.chi); end
    if ~exist('fill_value', 'var'), fill_value = nan; end

    if relation == '>'
        mask = maskvar > criterion;
    elseif relation == '<'
        mask = maskvar < criterion;
    elseif relation == '='
        mask = maskvar == criterion;
    end

    if isnan(fill_value)
        numnans = sum(isnan(chi.chi(Trange)));
    else
        numnans = sum(chi.chi(Trange) == fill_value);
    end

    chi.chi(mask) = fill_value;
    chi.eps(mask) = fill_value;
    chi.Kt(mask) = fill_value;
    chi.Jq(mask) = fill_value;
    chi.mask = isnan(chi.chi);

    if isnan(fill_value)
        numnewnans = sum(isnan(chi.chi(Trange)));
        verb = ' NaN-ed out ';
    else
        numnewnans = sum(chi.chi(Trange) == fill_value);
        if fill_value == 0
            verb = ' zeroed out ';
        else
            verb = [' filled (with ' num2str(fill_value) ') '];
        end
    end

    dnans = (numnewnans-numnans);
    percentage = dnans / length(Trange)*100;
    disp(['    ' maskname ' ' relation ' ' ...
          num2str(nanmedian(criterion), '%1.1e') ...
          verb num2str(percentage, '%.2f') ...
          '% of estimates (' num2str(dnans) ' points)'])

end
