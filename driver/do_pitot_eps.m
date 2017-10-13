%% this script does all the basic calibration excluding chi-processing and save 1s infomation
%     in ./proc/temp.mat
%  
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:01:26 PDT 2017


do_parallel = 1;     % use paralelle computing 
do_raw_proc = 1;     % do the raw data processing 
do_plot     = 0;     % generate a over view plot 

time_range = [datenum(2000, 1, 1, 0, 0, 0) ...
              datenum(2060, 1, 1, 0, 0, 0)];

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location


if do_raw_proc
   %_____________________get list of all raw data______________________
      [fids, fdate] = chi_find_rawfiles(basedir);


   %_____________processing loop through all raw files__________________

   disp('Calculating epsilon based on pitot data ')

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               proc_pitot_eps(basedir, fids{f});
               proc_pitot_eps(basedir, fids{f}, 2/(24*3600), [2.5 5]);
            catch ME
               disp(['!!!!!! ' fids{f} ' crashed while processing T structure !!!!!!' ]);
               disp(ME)
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            proc_pitot_eps(basedir, fids{f});
            proc_pitot_eps(basedir, fids{f}, 2/(24*3600), [2.5 5]);
         end
      end

   %____________________merge individual files______________________
      chi_merge_and_avg(basedir, 'pitot_eps_600sec', 0);;
      chi_merge_and_avg(basedir, 'pitot_eps_2sec', 0);
end


if do_plot
   load([basedir filesep 'proc' filesep 'pitot_eps.mat'])

	 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
	         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
	 
            [ax, ~] = create_axes(fig, 2, 1, 0);
            
            a=1;
            plot(ax(a), Peps.time, Peps.eps, 'Linewidth', 1);
            Peps.eps(abs(Peps.vel_avg)<.05) = nan;
            plot(ax(a), Peps.time, Peps.eps, 'Linewidth', 1);
            plot(ax(a), Peps.time, (2e-3*abs(Peps.vel_avg).^2).^1.5/.4, 'Linewidth', 1);
            plot(ax(a), Peps.time, Peps.eps_hf, 'Linewidth', 1);
            plot(ax(a), Peps.time, Peps.eps.*abs(Peps.vel_avg),'k', 'Linewidth', 1);


            legend('\epsilon_{pitot}', '\epsilon_{pitot>5cm/s}','\epsilon_{bbl1m}', '\epsilon_{pitotHF}',...
                    '\epsilon_{pitot}*|u|' );
            yl = 10.^[-9 -2];
            ylim(ax(a), yl);
            

            set(ax(a), 'Yscale', 'log');

            a=2;
            iif = [1:100 100:10:(length(Peps.f(:,1))*.5)] ;
            eps_ff = movmean(Peps.eps_f,3);
            pcolor(ax(a), Peps.time, Peps.f(iif,1), log10(eps_ff(iif,:)))
               shading(ax(a), 'flat');
               caxis(ax(a), [-8 -3]);
               cb = colorbar('peer', ax(a));
               cb.Position = [.95 .2 .01 .2];
               set(ax(a), 'Yscale', 'log');

               datetick(ax(a), 'keeplimits');
               
                     
               linkaxes(ax, 'x');
               
            
            print(gcf,'../pics/Pitot_eps.png','-dpng','-r200','-painters')
            
            

end
