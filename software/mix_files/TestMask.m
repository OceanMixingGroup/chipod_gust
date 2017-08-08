function [] = TestMask(chi, mask_var, mask_sign, mask_vals, name, iiTrange)

    hfig = figure;
    chiold = chi;
    Histograms(chi, hfig, 'count', 'raw');
    for ii=1:length(mask_vals)
        chi = ApplyMask(chi, mask_var, mask_sign, mask_vals(ii), ...
                        name, iiTrange);
        if mask_sign == '>'
            label = [name ' < ' num2str(mask_vals(ii), '%.1e')];
        else
            label = [name ' > ' num2str(mask_vals(ii), '%.1e')];
        end
        Histograms(chi, hfig, 'count', label);

        chi = chiold; % reset
    end