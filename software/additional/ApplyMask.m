function [chi, percentage, mask] = ApplyMask(chi, maskvar, relation, criterion, ...
                                             maskname, Trange, fill_value, ...
                                             do_plot, hfig, ID, label)

    if ~exist('Trange', 'var') | isempty(Trange), Trange = 1:length(chi.chi); end
    if ~exist('fill_value', 'var'), fill_value = nan; end
    if ~exist('do_plot', 'var'), do_plot = 0; end
    if do_plot
        if ~exist('hfig', 'var'), hfig = CreateFigure(); end
        if ~exist('ID', 'var'), ID=''; end
        if ~exist('label', 'var'), label = 'mask'; end
    end

    % make mask variable
    if relation == '>'
        mask = maskvar > criterion;
    elseif relation == '<'
        mask = maskvar < criterion;
    elseif relation == '='
        mask = maskvar == criterion;
    end

    if sum(mask) > 0 % fast path if mask is not going to do anything
        % find number of NaNs prior to masking
        if isnan(fill_value)
            numnans = sum(isnan(chi.chi(Trange)));
        else
            numnans = sum(chi.chi(Trange) == fill_value);
        end

        % apply the mask
        chi.chi(mask) = fill_value;
        if fill_value == 0
            warn_eps_nan = 1;
            chi.eps(mask) = NaN;
        else
            warn_eps_nan = 0;
            chi.eps(mask) = NaN;
        end
        chi.Kt(mask) = fill_value;
        chi.Jq(mask) = fill_value;
        chi.mask = isnan(chi.chi);

        % count number of NaNs (or fill_value) after masking
        if isnan(fill_value)
            numnewnans = sum(isnan(chi.chi(Trange)));
        else
            numnewnans = sum(chi.chi(Trange) == fill_value);
        end

        % save percentage and print to screen.
        dnans = (numnewnans-numnans);
        percentage = dnans / length(Trange)*100;

        if percentage > 0.5 & do_plot
            perlabel = [label ' -' num2str(percentage, '%.1f') '%'];
            Histograms(chi, hfig, 'count', ID, [perlabel]);
        end
    else
        percentage = 0;
        dnans = 0;
        warn_eps_nan = 0;
    end

    % get nice verb
    if isnan(fill_value)
        verb = ' NaN-ed out ';
    else
        if fill_value == 0
            verb = ' zeroed out ';
        else
            verb = [' filled (with ' num2str(fill_value) ') '];
        end
    end

    disp(['    ' maskname ' ' relation ' ' ...
          num2str(nanmedian(criterion), '%1.1e') ...
          verb num2str(percentage, '%.2f') ...
          '% of estimates (' num2str(dnans) ' points)'])

    if warn_eps_nan, disp('        Filling eps with NaN instead of 0.'); end

end
