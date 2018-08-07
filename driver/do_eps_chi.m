%%  This script runs the chi-processing based on Pitot epsilon in the BBL
%
%  
%   created by: 
%        Johannes Becherer
%        Fri Aug  3 11:53:11 PDT 2018


do_parallel = 0;     % use paralelle computing 
do_raw_proc = 0;     % do the raw data processing 
do_bbl_proc = 0;     % do processing basedon BBL scaling
    mab  = .5;        % height above the bottom
do_plot     = 1;     % generate a over view plot 


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    =   chi_get_unit_name(basedir); % get unit name
   rawdir  =   [basedir filesep 'raw' filesep]; % raw files location

% get time limits from whoAmI;
   [TL] =   whoAmI_timeLimits(basedir);
   time_range      = TL.master;

 % processing based on Pitot epsilon
if do_raw_proc

    load([basedir '/proc/pitot_eps.mat'])
    S.time  =  Peps.time;
    S.spd   =  Peps.spd;

    Eps.time   =  Peps.time;
    Eps.eps    =  Peps.eps;

    generate_eps_chi(basedir, Eps, S,  do_parallel,  time_range);


end

 %processing based on BBL scalinig
if do_bbl_proc   
    load([basedir '/proc/pitot_eps.mat'])
    S.time  =  Peps.time;
    S.spd   =  Peps.spd_lp;
    Eps.time   =  Peps.time;
    Eps.eps    =  (2e-3*abs(Peps.spd_lp).^2).^1.5/.4/mab;
    generate_eps_chi(basedir, Eps, S,  do_parallel,  time_range, '_bbl');
end



if do_plot

      ww = 120;


      load([basedir '/proc/eps_chi_bbl.mat']);

         chi.eps_m   =  movmean_nan( chi.eps, ww, .9);
         chi.chi_m   =  movmean_nan( chi.chi, ww, .9);
       [chi.Tz_eff, chi.Gamma, chi.Kt, chi.N2, chi.Loz]  =  Tz_bbl_eff( chi.eps, chi.chi, mab, chi.T, chi.S, chi.depth);
       [chi.Tz_eff_m, chi.Gamma_m, chi.Kt_m, chi.N2_m, chi.Loz_m]  =  Tz_bbl_eff( chi.eps_m, chi.chi_m, mab, chi.T, chi.S, chi.depth);
      chi_bbl   = chi;
         chi_bbl.chi(chi_bbl.chi<1e-14) = nan;
         chi_bbl.eps(chi_bbl.eps<1e-12) = nan;

      
      load([basedir '/proc/eps_chi.mat']);
         chi.eps_m   =  movmean_nan( chi.eps, ww, .9);
         chi.chi_m   =  movmean_nan( chi.chi, ww, .9);
       [chi.Tz_eff, chi.Gamma, chi.Kt, chi.N2, chi.Loz]  =  Tz_bbl_eff( chi.eps, chi.chi, mab, chi.T, chi.S, chi.depth);
       [chi.Tz_eff_m, chi.Gamma_m, chi.Kt_m, chi.N2_m, chi.Loz_m]  =  Tz_bbl_eff( chi.eps_m, chi.chi_m, mab, chi.T, chi.S, chi.depth);

      load([basedir '/input/dTdz_m.mat']);


   %_____________________Time series plots______________________
       fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
               'Papersize',[30 30],'PaperPosition',[0 0 30 30])
       [ax, ~] = create_axes(fig, 6, 1, 0);
       col = get(groot,'DefaultAxesColorOrder');
       clear p;
       
       tl = chi.time([1 end]);
       
       a=1;
       plot(ax(a), chi.time, chi.spd.^2,'color', [col(4,:) .5], 'Linewidth', 1);
       plot(ax(a), chi.time, movmean(chi.spd,ww).^2,'color', [col(4,:)*.5 1], 'Linewidth', 2);
        % soli mask
         tmp_msp2 =  (movmean(chi.spd,120).^2) - (movmean(chi.spd, 60*15).^2);
         %soli_mask   =  (tmp_msp2./movmean(chi.spd,120).^2)>.5 & tmp_msp2>.01;
         soli_mask   =  (tmp_msp2./movmean(chi.spd,60*20).^2)>1.5 & tmp_msp2>.01;
         soli_mask   =  movmean(soli_mask,180);
         plot(ax(a), chi.time, soli_mask*.1,'m', 'Linewidth', 2);
        % bore mask
         bore_mask   =  (movmean(chi.spd,30*60).^2)>.015;
         bore_mask   =  movmean(bore_mask,30*60);
         plot(ax(a), chi.time, bore_mask*.1,'c', 'Linewidth', 2);
         
         bore_head   = [(diff(bore_mask)>0) 0] == 1;
         plot(ax(a), chi.time, bore_head*.1,'b', 'Linewidth', 2);

         

            plot(ax(a), tl,[1 1]*0, 'k', 'Linewidth', 1);
         ylabel(ax(a), 'm^2/s^2')
         %yl = [-.1 .5];
         yl = [0 .2];
         ylim(ax(a), yl);
         xlim(ax(a), tl);
         t = text_corner(ax(a), ['u^2'], 1);
         t.FontWeight   =  'bold';
         t.BackgroundColor =  [1 1 1];

       a=2;
       plot(ax(a), chi_bbl.time, chi_bbl.chi,'color', [col(1,:) .5], 'Linewidth', 1);
       plot(ax(a), chi.time, chi.chi,'color', [col(2,:) .5], 'Linewidth', 1);
       p(1) = plot(ax(a), chi.time, chi.chi_m,'color', [col(2,:)*.5 1], 'Linewidth', 2);
       p(2) = plot(ax(a), chi_bbl.time, movmean_nan(chi_bbl.chi, ww, .9),'color', [col(1,:)*.5 1], 'Linewidth', 2);
       set(ax(a), 'Yscale', 'log');
         xlim(ax(a), tl);
         ylabel(ax(a), 'K^2/s')
         yl = [1e-12 1e-3];
         ylim(ax(a), yl);
       
			legend(p, 'pitot', 'bbl-scaling')
         t = text_corner(ax(a), ['\chi'], 1);
         t.FontWeight   =  'bold';
         t.BackgroundColor =  [1 1 1];
         
         
	

       a=3;
       plot(ax(a), chi.time, chi.eps,'color', [col(2,:) .5], 'Linewidth', 1);
       plot(ax(a), chi_bbl.time, chi_bbl.eps,'color', [col(1,:) .5], 'Linewidth', 2);
       p(1) = plot(ax(a), chi.time, chi.eps_m,'color', [col(2,:)*.5 1], 'Linewidth', 2);
       p(2) = plot(ax(a), chi_bbl.time, chi_bbl.eps_m,'color', [col(1,:)*.5 1], 'Linewidth', 2);
       set(ax(a), 'Yscale', 'log');
         xlim(ax(a), tl);
         ylabel(ax(a), 'm^2/s^3')
         yl = [1e-10 1e-3];
         ylim(ax(a), yl);
			legend(p, 'pitot', 'bbl-scaling')
         t = text_corner(ax(a), ['\epsilon'], 1);
         t.FontWeight   =  'bold';
         t.BackgroundColor =  [1 1 1];

       a=4;
         plot(ax(a), chi_bbl.time, chi_bbl.Tz_eff,'color', [col(1,:) .5], 'Linewidth', 1);
         plot(ax(a), chi.time, chi.Tz_eff,'color', [col(2,:) .5], 'Linewidth', 1);
         p(2) = plot(ax(a), chi_bbl.time, chi_bbl.Tz_eff_m,'color', [col(1,:)*.5 1], 'Linewidth', 2);
         p(1) = plot(ax(a), chi.time, chi.Tz_eff_m,'color', [col(2,:)*.5 1], 'Linewidth', 2);
         %plot(ax(a), chi.time, movmean_nan(chi.Tz_eff, ww, .9),'color',[.5 .5 .5 .5], 'Linewidth', 2);
         p(3)  =  plot(ax(a), Tz_m.time, Tz_m.Tz, 'color',[.3 .3 .3 .5], 'Linewidth', 2);
			legend(p, 'pitot', 'bbl-scaling', 'mooring')
         t = text_corner(ax(a), ['Tz_{eff}'], 1);
         t.FontWeight   =  'bold';
         t.BackgroundColor =  [1 1 1];
         yl = [1e-4 1e0];
         ylim(ax(a), yl);
         ylabel(ax(a), 'K/m')
         
         set(ax(a), 'Yscale', 'log');
         xlim(ax(a), tl);
         
       a=5;
         plot(ax(a), chi_bbl.time, chi_bbl.Gamma,'color', [col(1,:) .5], 'Linewidth', 1);
         plot(ax(a), chi.time, chi.Gamma,'color', [col(2,:) .5], 'Linewidth', 1);
         p(1) = plot(ax(a), chi.time, chi.Gamma_m,'color', [col(2,:)*.5 1], 'Linewidth', 2);
         p(2) = plot(ax(a), chi_bbl.time, chi_bbl.Gamma_m,'color', [col(1,:)*.5 1], 'Linewidth', 2);
         p(3) = plot(ax(a), chi.time, movmean_nan(chi.Gamma, ww, .9),'color', [.5 .5 .5 .5], 'Linewidth', 2);
            plot(ax(a), tl,[1 1]*.2, 'k', 'Linewidth', 1);
			legend(p, 'pitot', 'bbl-scaling', 'after mean pitot')
         yl = [1e-3 1e0];
         ylim(ax(a), yl);
         set(ax(a), 'Yscale', 'log');
         t = text_corner(ax(a), ['\Gamma'], 1);
         t.FontWeight   =  'bold';
         t.BackgroundColor =  [1 1 1];


       a=6;
         plot(ax(a), chi_bbl.time, chi_bbl.Loz/.4/mab,'color', [col(1,:) .5], 'Linewidth', 1);
         plot(ax(a), chi.time, chi.Loz/.4/mab,'color', [col(2,:) .5], 'Linewidth', 1);
         p(1) = plot(ax(a), chi.time, chi.Loz_m/.4/mab,'color', [col(2,:)*.5 1], 'Linewidth', 2);
         p(2) = plot(ax(a), chi_bbl.time, chi_bbl.Loz_m/.4/mab,'color', [col(1,:)*.5 1], 'Linewidth', 2);
			   legend(p, 'pitot', 'bbl-scaling')
            plot(ax(a), tl,[1 1]*5^(3/4), 'k', 'Linewidth', 1);
            yl = [1e0 1e3];
            ylim(ax(a), yl);
            set(ax(a), 'Yscale', 'log');
            t = text_corner(ax(a), ['L_{Oz}/({\kappa}z)'], 1);
            t.FontWeight   =  'bold';
            t.BackgroundColor =  [1 1 1];

       linkaxes(ax, 'x');
         xlim(ax(a), tl);
       datetick(ax(a), 'keeplimits');

       print(gcf,'../pics/eps_fit_chi.png','-dpng','-r200','-painters')
       
       
%%_____________________correlation 2d hist______________________
    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[30 12],'PaperPosition',[0 0 30 12])
    
       load cmap;
            [ax, ~] = create_axes(fig, 1, 4, 0);
               squeeze_axes(ax, 1,.95)


        % epsilon    
            a=1;
             chi.eps_bbl = chi_bbl.eps_m;
             chi.eps_m( chi.eps_m <=0) = nan;
             chi.eps_bbl(isnan(chi.eps_m)) =  nan;

             X   = log10(chi_bbl.eps_m);
             Y   = log10(chi.eps_m);
              relativeQ =   sum(~isnan(chi.eps_m))./length(chi.eps_m)*100;
               XL = 'BBL';
               YL = 'pitot';
               xl = [-7.5 -4.5];
               binx = [xl(1):(diff(xl)/20):xl(2)];

               plot_hist2d(ax(a), X, Y, binx, binx,  XL, YL, 1 )

               t = text_corner(ax(a), ['Good Data (' num2str(round(relativeQ)) '% of all data)'], -2);
               t.FontWeight  =  'bold';
               t = text_corner(ax(a), ['log_{10}\epsilon'], 1);
               t.FontSize  = 18;

           a=2;
            X   = log10(chi_bbl.eps_m(bore_mask>0));
            Y   = log10(chi.eps_m(bore_mask>0));
              relativeQ =   sum(bore_mask>0 & ~isnan(chi.eps_m))./sum(~isnan(chi.eps_m))*100;
            
            plot_hist2d(ax(a), X, Y, binx, binx,  XL, YL, 1 )
               t = text_corner(ax(a), ['bore (' num2str(round(relativeQ)) '% of good data)'], -2);
               t.FontWeight  =  'bold';
               t = text_corner(ax(a), ['log_{10}\epsilon'], 1);
               t.FontSize  = 18;
            
       
           a=3;
            X   = log10(chi_bbl.eps_m(soli_mask>0));
            Y   = log10(chi.eps_m(soli_mask>0));
              relativeQ =   sum(soli_mask>0 & ~isnan(chi.eps_m))./sum(~isnan(chi.eps_m))*100;
            plot_hist2d(ax(a), X, Y, binx, binx,  XL, YL, 1 )
               t = text_corner(ax(a), ['soliton (' num2str(round(relativeQ)) '% of good data)'], -2);
               t.FontWeight  =  'bold';
               t = text_corner(ax(a), ['log_{10}\epsilon'], 1);
               t.FontSize  = 18;
           a=4;
            X   = log10(chi_bbl.eps_m(bore_head>0));
            Y   = log10(chi.eps_m(bore_head>0));
              relativeQ =   sum(bore_head>0 & ~isnan(chi.eps_m))./sum(~isnan(chi.eps_m))*100;
            plot_hist2d(ax(a), X, Y, binx, binx,  XL, YL, 1 )
               t = text_corner(ax(a), ['bore_head (' num2str(round(relativeQ)) '% of good data)'], -2);
               t.FontWeight  =  'bold';
               t = text_corner(ax(a), ['log_{10}\epsilon'], 1);
               t.FontSize  = 18;
       print(gcf,'../pics/BBL_vs_pitot_bore_soliton.png','-dpng','-r200','-painters')
          
%%_____________________correlation 2d hist______________________
    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[30 10],'PaperPosition',[0 0 30 10])
    
       load cmap;
            [ax, ~] = create_axes(fig, 1, 3, 0);
               squeeze_axes(ax, 1,.9)
               shift_axes(ax, 0, .05)
            clear X;
         a=1;
            X{1}   = log10(chi.eps_m);
            X{2}   = log10(chi.eps_m(bore_mask>0));
            X{3}   = log10(chi.eps_m(soli_mask>0));
            sl = [min(X{1}) max(X{1})];
            bins = [sl(1):(diff(sl)/50):sl(2)];
         p = plot_hist1d(ax(a), X, bins, 'log_{10}\epsilon', 0)
         legend(p, 'All', 'bore', 'soliton');
         
         a=2;
             chi.chi_m( chi.chi_m <=0) = nan;
            X{1}   = log10(chi.chi_m);
            X{2}   = log10(chi.chi_m(bore_mask>0));
            X{3}   = log10(chi.chi_m(soli_mask>0));
            sl = [min(X{1}) max(X{1})];
            bins = [sl(1):(diff(sl)/50):sl(2)];
         p = plot_hist1d(ax(a), X, bins, 'log_{10}\chi', 0)
         legend(p, 'All', 'bore', 'soliton');

         a=3;
            X{1}   = log10(chi.Gamma_m);
            X{2}   = log10(chi.Gamma_m(bore_mask>0));
            X{3}   = log10(chi.Gamma_m(soli_mask>0));
            X{4}   = log10(chi.Gamma_m(bore_head>0));
            sl = [min(X{1}) max(X{1})];
            bins = [sl(1):(diff(sl)/50):sl(2)];
         [p, Ncnt] = plot_hist1d(ax(a), X, bins, 'log_{10}\Gamma', 0)
            yl = [0 1.05]*max([Ncnt{:}]);
            plot(ax(a), [1 1]*log10(.2), yl, 'k--',  'Linewidth', 1);
            
         legend(p, 'All', 'bore', 'soliton', 'bore-head');


       print(gcf,'../pics/bore_soliton.png','-dpng','-r200','-painters')

    % [fig] = plot_pitot_eps( basedir)
    % print(fig,['../pics/Pitot_eps.png'],'-dpng','-r200','-painters')
end

function [p, Ncnt] = plot_hist1d(ax, X, bins, label, normalized)
   col = get(groot,'DefaultAxesColorOrder');
      sl = bins([1 end]);
      for i =1:length(X);
         [Ncnt{i},~] = histcounts( X{i} , bins);
         if normalized
            Ncnt{i} = Ncnt{i}/max([Ncnt{i}]);
         end
      end
      yl = [0 1.05]*max([Ncnt{:}]);
      for i =1:length(X);
         pj = i; p(pj) = plot( ax, bins(1:end-1)+diff(bins(1:2)*.5), ...
                              Ncnt{i} , 'color', [col(pj,:) 1], 'Linewidth', 2);   
         patch( [bins(1) bins(1:end)]+diff(bins(1:2)*.5), [0 Ncnt{i} 0], ...
                  col(i,:)*.5+.5,'Facealpha',.3, 'Parent', ax);
         plot(ax, nanmedian(X{i}), yl(2), '+', 'color', [col(pj,:) 1]);
         plot(ax, nanmean(X{i}), yl(2), '*', 'color', [col(pj,:) 1]);
      end
      xlabel(ax, [label])
      %t = text_corner(ax, ['\langle label \rangle = ' num2str(nanmean(X{1}), '%3.1e') ], 36);
      %t.Color = col(1,:);
      set(ax, 'Ycolor', [1 1 1]);
      xlim(ax, sl);
      ylim(ax, yl);
end

function [ax] = plot_hist2d(ax, X, Y, binx, biny,  XL, YL, gridOn )
%% [ax] = plot_hist2d(ax, X, Y, binx, biny,  XL, YL, gridOn )
%        
%     This function plots a 2 d histogram in the ax
%
%     INPUT
%        ax    :  axis handle
%        X     :  signal 1;
%        Y     :  signal 2:
%        binx  :  bins for 1
%        biny  :  bins for 2
%        XL    :  label for X
%        YL    :  label for Y
%        gridOn:  should I plot a grid ('1' default)
%
%
%   created by: 
%        Johannes Becherer
%        Mon Aug  6 11:50:46 PDT 2018

xl = binx([1 end]);
yl = biny([1 end]);

      load cmap;
            [hist,~,~,~] = hist2d(binx, biny, X, 0, Y, 0, 3);
            %pcolor(ax, binx, biny, hist);    shading(ax,'flat');
            contourf(ax,binx,biny,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
               if gridOn
                  plot(ax, binx, biny,'k', 'Linewidth', 1);
                  plot(ax, binx, log10(5*10.^biny),'k:', 'Linewidth', 1);
                  plot(ax, binx, log10(2*10.^biny),'k--', 'Linewidth', 1);
                  plot(ax, binx, log10(.5*10.^biny),'k--', 'Linewidth', 1);
                  plot(ax, binx, log10(.2*10.^biny),'k:', 'Linewidth', 1);
               end
               colormap(ax, cmap.chi);
               t = text_corner(ax, [XL], 8);
                  t.BackgroundColor = [1 1 1];
                  t.FontWeight      = 'bold';
               t = text_corner(ax, [YL], 4);
                  t.BackgroundColor = [1 1 1];
                  t.FontWeight      = 'bold';
               ylim(ax, yl);
               xlim(ax, xl);
               % corelation
                  if size(X,2)>size(X,1)
                     X = X';
                  end
                  if size(Y,2)>size(Y,1)
                     Y = Y';
                  end
                  [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
                  t = text_corner(ax,  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                           ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);
               set(ax, 'box', 'on')
end


function [x_mean]   = movmean_nan( x, k, maxPnan)
      x_mean   =  movmean( x, k, 'omitnan').*(movmean(isnan(x),k)<maxPnan);
end
