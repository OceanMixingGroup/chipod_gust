function [fig_noise_floor] = make_noise_floor_histograms(...
    ID, chi, fig_noise_floor, is_visible, spec_floor, spec_floor_mask, do_spec_area)

    if isempty(fig_noise_floor)
        fig_noise_floor.fig = CreateFigure(is_visible, ...
                                           'Noise floor diagnostics');
        fig_noise_floor.fig.Position(3) = 2000;
        fig_noise_floor.fig.Position(4) = 1065;
        fig_noise_floor.tparea = subplot(2,2,[1,2]); hold on;
        fig_noise_floor.chi = subplot(223); hold on;
        fig_noise_floor.eps = subplot(224); hold on;
    end

    clr = choose_color(ID, 'color');

    if do_spec_area
        % histogram of area under sepctrum
        histogram(fig_noise_floor.tparea, log10(chi.spec_area), 'EdgeColor', clr, ...
                  'displaystyle', 'stairs', 'DisplayName', ID);
        % area under spectrum if spectrum was noise
        plot(fig_noise_floor.tparea, ...
             [1 1]* log10(nanmedian(spec_floor) * nanmean(chi.nfft)), ...
             fig_noise_floor.tparea.YLim, ...
             '--', 'color', clr, 'displayname', [ID ' noise floor'])
        legend(fig_noise_floor.tparea, '-dynamiclegend');
        xlabel(fig_noise_floor.tparea, 'log_{10} area under Tp spectrum')
        ylabel(fig_noise_floor.tparea, 'count')
    end

    histogram(fig_noise_floor.chi, log10(chi.chi(~spec_floor_mask)), 'EdgeColor', clr, ...
              'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(fig_noise_floor.chi, log10(chi.chi(spec_floor_mask)), 'linestyle', '--', ...
              'EdgeColor', clr, 'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(fig_noise_floor.chi, log10(chi.chi), 'EdgeColor', clr, 'linewidth', 2, ...
              'normalization', 'count', 'displaystyle', 'stairs', 'displayname', ['\chi_{' ID '}'])
    legend(fig_noise_floor.chi, '-dynamiclegend')
    xlabel(fig_noise_floor.chi, 'log_{10} \chi')
    ylabel(fig_noise_floor.chi, 'count')
    title(fig_noise_floor.chi, ...
          {'Fits can work with white noise input.';
           'thin dashed, solid = successful fits (below, above) noise floor';
           'thick = all values'})

    histogram(fig_noise_floor.eps, log10(chi.eps(~spec_floor_mask)), 'EdgeColor', clr, ...
              'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(fig_noise_floor.eps, log10(chi.eps(spec_floor_mask)), 'linestyle', '--', ...
              'EdgeColor', clr, 'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(fig_noise_floor.eps, log10(chi.eps), 'EdgeColor', clr, 'linewidth', 2, ...
              'normalization', 'count', 'displaystyle', 'stairs', 'displayname', ['\epsilon_{' ID '}'])
    legend(fig_noise_floor.eps, '-dynamiclegend')
    xlabel(fig_noise_floor.eps, 'log_{10} \epsilon')
    ylabel(fig_noise_floor.eps, 'count')
    title(fig_noise_floor.eps, ...
          {'Fits can work with white noise input.';
           'thin dashed, solid = successful fits (below, above) noise floor';
           'thick = all values'})

end
