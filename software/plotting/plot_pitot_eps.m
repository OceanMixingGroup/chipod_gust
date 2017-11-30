function [fig] = plot_pitot_eps( basedir, spec_length)
%%   [fig] = plot_pitot_eps( basedir, [spec_length])
%
%  This function generate a basic diagnostic plot for the epsilon from Pitot tube
%
%  INPUT
%     basedir   :  instrument directory
%     matfile   :  which mat file (default '300sec')
%
%   created by: 
%        Johannes Becherer
%        Thu Nov 30 14:54:08 PST 2017

if nargin < 2
   spec_length = '300sec';
end

   unit    = chi_get_unit_name(basedir); % get unit name

   load([basedir filesep 'proc' filesep 'pitot_eps' spec_length  '.mat'])

	 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
	         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
	 
            [ax, ~] = create_axes(fig, 2, 1, 0);
            col = get(groot,'DefaultAxesColorOrder');
            
            xl = Peps.time([1 end]);
            
            Navg = round(length(Peps.time)/100);
            if Navg<1
               Navg=1;
            end
            
            a=1;
               plot(ax(a), Peps.time, Peps.eps, 'color', [col(1,:) .5], 'Linewidth', 1);
               pj = 1; p(pj) = plot(ax(a), Peps.time, movmean(Peps.eps, Navg, 'omitnan'), 'color', [col(pj,:)*.5 1], 'Linewidth', 2);
               plot(ax(a), Peps.time, (2e-3*abs(Peps.vel).^2).^1.5/.4, 'color', [col(2,:) .5], 'Linewidth', 1);
               pj = 2; p(pj) =plot(ax(a), Peps.time, movmean((2e-3*abs(Peps.vel).^2).^1.5/.4, Navg, 'omitnan'), 'color', [col(2,:)*.7 1], 'Linewidth', 2);
               

               legend(p, '\epsilon_{pitot}', '\epsilon_{bbl1m}');
               yl = 10.^[-9 -2];
               ylim(ax(a), yl);
               xlim(ax(a), xl);
               t = text_corner(ax(a), ['\epsilon [m2/s3]'], 1);
               t = text_corner(ax(a), [unit], -2);
               
               
               
               
               
               set(ax(a), 'Yscale', 'log');

            a=2;
               plot(ax(a), Peps.time, Peps.vel, 'color', [col(1,:) .5], 'Linewidth', 1);
               plot(ax(a), Peps.time,  movmean(Peps.vel, Navg, 'omitnan'), 'color', [col(1,:)*.5 1], 'Linewidth', 2);
               t = text_corner(ax(a), ['speed [m/s]'], 1);

               xlim(ax(a), xl);
               datetick(ax(a), 'keeplimits');
               
                     
               linkaxes(ax, 'x');


return
   % spectro gram
     %      iif = [1:100 100:10:(length(Peps.f(:,1))*.5)] ;
     %      eps_ff = movmean(Peps.eps_f,3);
     %      pcolor(ax(a), Peps.time, Peps.f(iif,1), log10(eps_ff(iif,:)))
     %         shading(ax(a), 'flat');
     %         caxis(ax(a), [-8 -3]);
     %         cb = colorbar('peer', ax(a));
     %         cb.Position = [.95 .2 .01 .2];
     %         set(ax(a), 'Yscale', 'log');

     %         
     %      
     %      print(gcf,'../pics/Pitot_eps.png','-dpng','-r200','-painters')
            
            

