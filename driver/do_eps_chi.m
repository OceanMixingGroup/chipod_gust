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
       [chi.Tz_eff, chi.Gamma, chi.Kt]  =  Tz_bbl_eff( chi.eps, chi.chi, mab, chi.T, chi.S, chi.depth);
       [chi.Tz_eff_m, chi.Gamma_m, chi.Kt_m]  =  Tz_bbl_eff( chi.eps_m, chi.chi_m, mab, chi.T, chi.S, chi.depth);
      chi_bbl   = chi;
         chi_bbl.chi(chi_bbl.chi<1e-14) = nan;
         chi_bbl.eps(chi_bbl.eps<1e-12) = nan;

      
      load([basedir '/proc/eps_chi.mat']);
         chi.eps_m   =  movmean_nan( chi.eps, ww, .9);
         chi.chi_m   =  movmean_nan( chi.chi, ww, .9);
       [chi.Tz_eff, chi.Gamma, chi.Kt]  =  Tz_bbl_eff( chi.eps, chi.chi, mab, chi.T, chi.S, chi.depth);
       [chi.Tz_eff_m, chi.Gamma_m, chi.Kt_m]  =  Tz_bbl_eff( chi.eps_m, chi.chi_m, mab, chi.T, chi.S, chi.depth);

      load([basedir '/input/dTdz_m.mat']);

       fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
               'Papersize',[30 20],'PaperPosition',[0 0 30 20])
       [ax, ~] = create_axes(fig, 5, 1, 0);
       col = get(groot,'DefaultAxesColorOrder');
       clear p;
       
       tl = chi.time([1 end]);
       
       a=1;
       plot(ax(a), chi.time, chi.spd.^2,'color', [col(4,:) .5], 'Linewidth', 1);
       plot(ax(a), chi.time, movmean(chi.spd,ww).^2,'color', [col(4,:)*.5 1], 'Linewidth', 2);
        % soli mask
         tmp_msp2 =  (movmean(chi.spd,120).^2) - (movmean(chi.spd, 60*30).^2);
         soli_mask   =  (tmp_msp2./movmean(chi.spd,120).^2)>.7 & tmp_msp2>.01;
         soli_mask   =  movmean(soli_mask,180);
         plot(ax(a), chi.time, soli_mask*.1,'m', 'Linewidth', 2);
        % bore mask
         bore_mask   =  (movmean(chi.spd,30*60).^2)>.02;
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

       linkaxes(ax, 'x');
         xlim(ax(a), tl);
       datetick(ax(a), 'keeplimits');

       
%%_____________________correlation______________________
    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[20 30],'PaperPosition',[0 0 20 30])
    
       load cmap;
            [ax, ~] = create_axes(fig, 1, 3, 0);
            
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
               yl = xl;
               binx = [xl(1):(diff(xl)/20):xl(2)];
               biny = binx;
            [hist,~,~,~] = hist2d(binx, biny, X, 0, Y, 0, 3);
            %pcolor(ax(a), binx, biny, hist);    shading(ax(a),'flat');
            contourf(ax(a),binx,biny,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
               plot(ax(a), binx, biny,'k', 'Linewidth', 1);
               colormap(ax(a), cmap.chi);
               xlabel(ax(a), XL)
               ylabel(ax(a), YL)
               ylim(ax(a), yl);
               xlim(ax(a), xl);
               % corelation
                  if size(X,2)>size(X,1)
                     X = X';
                  end
                  if size(Y,2)>size(Y,1)
                     Y = Y';
                  end
                  [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
                  t = text_corner(ax(a),  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                           ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);
               t = text_corner(ax(a), ['Good Data (' num2str(round(relativeQ)) '% of all data)'], 1);
               t.FontWeight  =  'bold';
                           
                           
            
           a=2;
            X   = log10(chi_bbl.eps_m(bore_mask>0));
            Y   = log10(chi.eps_m(bore_mask>0));
              relativeQ =   sum(bore_mask>0 & ~isnan(chi.eps_m))./sum(~isnan(chi.eps_m))*100;
            [hist,~,~,~] = hist2d(binx, biny, X, 0, Y, 0, 3);
            %pcolor(ax(a), binx, biny, hist);    shading(ax(a),'flat');
            contourf(ax(a),binx,biny,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
               colormap(ax(a), cmap.chi);
               plot(ax(a), binx, biny,'k', 'Linewidth', 1);
               xlabel(ax(a), XL)
               %ylabel(ax(a), YL)
               ylim(ax(a), yl);
               xlim(ax(a), xl);
                  [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
                  t = text_corner(ax(a),  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                           ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);
               t = text_corner(ax(a), ['bore (' num2str(round(relativeQ)) '% of good data)'], 1);
               t.FontWeight  =  'bold';
            
       
           a=3;
            X   = log10(chi_bbl.eps_m(soli_mask>0));
            Y   = log10(chi.eps_m(soli_mask>0));
              relativeQ =   sum(soli_mask>0 & ~isnan(chi.eps_m))./sum(~isnan(chi.eps_m))*100;
            [hist,~,~,~] = hist2d(binx, biny, X, 0, Y, 0, 3);
            %pcolor(ax(a), binx, biny, hist);    shading(ax(a),'flat');
            contourf(ax(a),binx,biny,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
               colormap(ax(a), cmap.chi);
               plot(ax(a), binx, biny,'k', 'Linewidth', 1);
               xlabel(ax(a), XL)
               %ylabel(ax(a), YL)
               ylim(ax(a), yl);
               xlim(ax(a), xl);
                  [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
                  t = text_corner(ax(a),  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                           ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);
               t = text_corner(ax(a), ['soliton (' num2str(round(relativeQ)) '% of good data)'], 1);
               t.FontWeight  =  'bold';
          


    % [fig] = plot_pitot_eps( basedir)
    % print(fig,['../pics/Pitot_eps.png'],'-dpng','-r200','-painters')
end

function [x_mean]   = movmean_nan( x, k, maxPnan)
      x_mean   =  movmean( x, k, 'omitnan').*(movmean(isnan(x),k)<maxPnan);
end
