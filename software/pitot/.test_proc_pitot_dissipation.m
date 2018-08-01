clear all;
close all;

addpath(genpath('./chipod_gust/software/'));

basedir = '/home/johannes/ismoor/oc50_a/gust/G031/';
basedir = '/home/johannes/ismoor/oc40s_t/gust/G019/';
%rfid    = 'raw_1709150000.G031';
rfid    = 'raw_1709230000.G031';
rfid    = 'raw_1709190000.G019';


proc_pitot_dissipation(basedir, rfid, 1/3600/24, [1 20], 1)


load ~/ismoor/oc40s_t/gust/G019/proc/pitot_dissipation/pitot_eps_1709190000.mat
%load ~/ismoor/oc50_a/gust/G031/proc/pitot_dissipation/pitot_eps_1709230000.mat
%load ~/ismoor/oc50_a/gust/G031/proc/pitot_dissipation/pitot_eps_1709150000.mat
   [Peps1]  = flag_Peps(Peps)
   Peps1.eps_bbl   =  (2e-3*Peps1.spd_lp.^2).^1.5/.4;

load ~/ismoor/oc40s_t/gust/G019/proc/chi/chi_mmg/chi_1709190000.mat
%load ~/ismoor/oc50_a/gust/G031/proc/pitot_eps2sec/pitot_eps_2sec_1709230000.mat 
%load ~/ismoor/oc50_a/gust/G031/proc/pitot_eps2sec/pitot_eps_2sec_1709150000.mat 
   %[Peps]  = flag_Peps(Peps)
   Peps  =  chi;
   %Peps.eps = Peps.eps*(2*pi)^2;
   Peps.eps_bbl   =  (2e-3*Peps.spd.^2).^1.5/.4;
   %Peps.eps_bbl   =  (2e-3*Peps.spd_lp.^2).^1.5/.4;



figure
   for a=1:2
      ax(a) = subplot( 2, 1, a);
      hold(ax(a), 'on');
   end
   a=1;
   plot(ax(a), Peps1.time, Peps1.eps, 'Linewidth', 1);
   plot(ax(a), Peps.time, Peps.eps, 'Linewidth', 1);
   plot(ax(a), Peps1.time, Peps1.eps_bbl, 'Linewidth', 1);
   plot(ax(a), Peps1.time, movmean(Peps1.eps,60, 'omitnan'), 'Linewidth', 1);
   set(ax(a), 'Yscale', 'log');

   a=2;
   plot(ax(a), Peps1.time, Peps1.spd, 'Linewidth', 1);
   plot(ax(a), Peps.time, Peps.spd, 'Linewidth', 1);
   linkaxes(ax, 'x');

figure
   subplot(1,2,1);
   loglog(movmean(Peps1.eps,60, 'omitnan'), Peps1.eps_bbl, '.');
   hold all;
   plot( [1e-10 1e-2],[1e-10 1e-2],'k',  'Linewidth', 1);

   a=1;
   col = get(groot,'DefaultAxesColorOrder');
   ax(1) = subplot(1,2,2)
   hold all;
   X{1}   = log10(abs(Peps1.eps));
   X{2}   = log10(abs(Peps1.eps_bbl));
   X{3}   = log10(abs(Peps.eps));
   sl = [min(X{1}) max(X{1})];
   bins = [sl(1):(diff(sl)/20):sl(2)];
   for i =1:length(X);
      [Ncnt{i},~] = histcounts( X{i} , bins);
      Ncnt{i} = Ncnt{i}/max(Ncnt{i});
      pj = i; p(pj) = plot(gca, bins(1:end-1)+diff(bins(1:2)*.5), ...
                           Ncnt{i} , 'color', [col(pj,:) 1], 'Linewidth', 2);   
      patch( [bins(1) bins(1:end)]+diff(bins(1:2)*.5), [0 Ncnt{i} 0], ...
               col(i,:)*.5+.5,'Facealpha',.3, 'Parent', ax(a));
   end
   xlabel(ax(a), ['label'])
   %t = text_corner(ax(a), ['\langle label \rangle = ' num2str(nanmean(X{1}), '%3.1e') ], 36);
   %t.Color = col(1,:);
   set(ax(a), 'Ycolor', [1 1 1]);
   xlim(ax(a), sl);
   
