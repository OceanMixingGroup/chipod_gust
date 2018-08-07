clear all;
close all;


load ../proc/pitot_eps

f = Peps.f_range(1):diff(Peps.f_range([1 2]))/(size(Peps.k,1)-1):Peps.f_range(2);

 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
 
            ax(1) = subplot( 3, 1, 1);
            hold(ax(1), 'on');
            ax(2) = subplot( 3, 1, [2 3]);
            hold(ax(2), 'on');

            
            ii = [1:1:10000]  + 6e4;

            a=1;
            plot(ax(a), Peps.time(ii), Peps.eps(ii));
               set(ax(a), 'Yscale', 'log');

            a=2;


            pcolor(ax(a), Peps.time(ii), f, log10(Peps.D_k(:,ii)));
               plot(ax(a), Peps.time(ii), Peps.spd(ii)./Peps.eta(ii), 'k')
               shading(ax(a), 'flat');
               colormap(jet);
               caxis([-11 -8]);
               set(ax(a), 'Yscale', 'log');
               yl = Peps.f_range([1 2]);
               ylim(ax(a), yl);
               

            linkaxes(ax, 'x');
            datetick(ax(a), 'keeplimits');
            
            

figure
   eps_target =  1e-6;
   ii6 = find( Peps.eps>.90*eps_target & Peps.eps<1.1*eps_target );
   
      k_test = reshape(Peps.k(:,ii6),1,numel(Peps.k(:,ii6)));
      D_test = reshape(Peps.D_k(:,ii6),1,numel(Peps.D_k(:,ii6)));
         [k_test ii_sort] = sort(k_test);
         D_test = D_test(ii_sort);

      plot( k_test, D_test);
      hold(gca,'on');
      plot( k_test, movmean(D_test,100), 'Linewidth',3);
      k_pao =  1e-1:1e-1:3e3;
      plot( k_pao, pao_spectrum(k_pao, eps_target));
      set(gca, 'Yscale', 'log', 'Xscale', 'log');
   
