clear all;
%close all;

addpath(genpath('./chipod_gust/software/'));

if 1 
   load ../proc/pitot_eps.mat;
end

Peps.eps = Peps.eps_nomask.*Peps.mask_tot;
nu = 1.6e-6;
eps(1) = 1e-7;
eps(2) = 3e-7;
eps(3) = 1e-6;
eps(4) = 3e-6;
eps(5) = 1e-5;
% nathsyth
   [G1,kks]=nasmyth_G1(1000,20);
   dk_grid = 6;
   k_grid      = 10:dk_grid:500;
   D_grid      = nan(length(k_grid),length(eps));
   D_grid_acc = nan(length(k_grid),length(eps));
   for i=1:length(eps)
      [D_na{i}, k_na{i}, eta(i)] = nasmyth2Dk( kks, G1, eps(i), nu );
      ii_eps = find(Peps.eps>10^(log10(eps(i))-.1) & Peps.eps<10^(log10(eps(i))+.1));
      tmp = [Peps.k(:,ii_eps)];
      k_vec = tmp(:);
      tmp = [Peps.D_k(:,ii_eps)];
      D_vec = tmp(:);
      tmp = [Peps.D_k_acc(:,ii_eps)];
      D_vec_acc = tmp(:);
      for j = 1:length(k_grid)
         ii_kbin   =  find( k_vec>=(k_grid(j)-.5*dk_grid) & k_vec<=(k_grid(j)+.5*dk_grid)  );
         D_grid(j,i)    = nanmean(D_vec(ii_kbin));
         D_grid_acc(j,i) = nanmean(D_vec_acc(ii_kbin));
         if length(ii_kbin)<10
            D_grid(j,i)= nan;
            D_grid_acc(j,i)= nan;
         end
      end
   end



 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
         'Papersize',[12 10],'PaperPosition',[0 0 12 10])
   
         for a=1:1
            ax(a) = subplot( 1, 1, a);
            hold(ax(a), 'on');
         end

         a=1;
         %plot(ax(a), k_vec, D_vec, '.', 'color', [.7 .7 .7 .5 ]);
         col = get(groot,'DefaultAxesColorOrder');
         
         xl = [6e-1 1e3];
         xlim(ax(a), xl);
         yl = [1e-6 3e-3];
         ylim(ax(a), yl);
         for i =1:length(eps)    
            plot(ax(a), k_grid, D_grid(:,i),'color',col(i,:)*.2+.7, 'Linewidth', 1);
            ii_kmax = find(k_grid<(1/(3*eta(i))));
            plot(ax(a), k_grid(ii_kmax), D_grid(ii_kmax,i),'color',col(i,:), 'Linewidth', 2);
            plot(ax(a), k_na{i}, D_na{i}, '--', 'color',col(i,:)*.5, 'Linewidth', 1);
            plot(ax(a), k_grid(ii_kmax), D_grid_acc(ii_kmax,i),'--','color',col(i,:), 'Linewidth', 1);

            i_text = min(find(k_na{i}>xl(1)));
            t= text( xl(1)*1.3, D_na{i}(i_text), ['\epsilon = ' num2str(eps(i)) ' m^2/s^3']);
            t.Color = col(i,:)*.8;
            %t.FontWeight = 'bold';
            t.BackgroundColor = [1 1 1 .8];
            
         end
         ylabel(ax(a), 'D(k) [s^{-2}(rad/m)^{-1}]');
         xlabel(ax(a), 'k [(rad/m)]');

         set(ax(a), 'Yscale', 'log', 'Xscale', 'log');
            set(ax(a), 'box', 'on', 'TickDir', 'out');

         print(gcf,'../pics/pitot_spec.png','-dpng','-r200','-painters')
         
         
