%% this script does all the basic calibration excluding chi-processing and save 1s infomation
%     in ./proc/temp.mat
%  
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:01:26 PDT 2017


do_parallel = 1;     % use paralelle computing 
do_raw_proc = 1;     % do the raw data processing 
do_plot     = 1;     % generate a over view plot 



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

   disp('calibrating all base quantities in ./proc/temp.mat ')

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               chi_T_proc(basedir, fids{f});
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
            chi_T_proc(basedir, fids{f});
         end
      end

   %____________________merge individual files______________________
      % average 1 sec
      chi_merge_and_avg(basedir, 'temp', 0);
end


%_____________________do an over view plot______________________
if do_plot
     load([basedir 'proc' filesep 'temp.mat']);

   tl = T.time([1 end]);


   col = get(groot,'DefaultAxesColorOrder');
   

    fig = figure('Color',[1 1 1],'visible',vis,'Paperunits','centimeters',...
            'Papersize',[30 30],'PaperPosition',[0 -1 30 30])


            [ax, ~] = create_axes(fig, 6, 1, 0);
            
            a = 1;
            pj = 1; p(pj) = plot(ax(a), T.time, T.depth, 'color', [0  0 0 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['P [dbar]'], 7);
               
            a = 2;
            if isfield(T, 'T1') 
               pj = 1; p(pj) = plot(ax(a), T.time, T.T1, 'color', [col(pj,:) 1], 'Linewidth', 2);
               pj = 2; p(pj) = plot(ax(a), T.time, T.T2, 'color', [col(pj,:) 1], 'Linewidth', 2);
                  t1 = text_corner(ax(a), ['T1 [^\circ C]'], 1);
                  t1.Color = [col(1,:)];
                  t2 = text_corner(ax(a), ['T2 [^\circ C]'], 3);
                  t2.Color = [col(2,:)];
            else
               pj = 1; p(pj) = plot(ax(a), T.time, T.T, 'color', [0  0 0 1], 'Linewidth', 2);
                  xlim(ax(a), tl);
                  t = text_corner(ax(a), ['T [^\circ C]'], 1);
            end
               xlim(ax(a), tl);
               
            
            a = 3;
            pj = 1; p(pj) = plot(ax(a), T.time, T.AX, 'color', [col(pj,:) 1], 'Linewidth', 2);
            pj = 2; p(pj) = plot(ax(a), T.time, T.AY, 'color', [col(pj,:) 1], 'Linewidth', 2);
            pj = 3; p(pj) = plot(ax(a), T.time, T.AZ, 'color', [col(pj,:) 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t1 = text_corner(ax(a), ['AX [m s^{-2}]'], 1);
               t1.Color = [col(1,:)];
               t2 = text_corner(ax(a), ['AY [m s^{-2}]'], 2);
               t2.Color = [col(2,:)];
               t3 = text_corner(ax(a), ['AZ [m s^{-2}]'], 3);
               t3.Color = [col(3,:)];
            
            a = 4;
            pj = 1; p(pj) = plot(ax(a), T.time, T.cmp, '.', 'color', [0  0 0 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['compass [deg]'], 7);

            a = 5;
            pj = 1; p(pj) = plot(ax(a), T.time, T.W, 'color', [0  0 0 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['pitot [volt]'], 7);

            a = 6;
            pj = 1; p(pj) = plot(ax(a), T.time, T.T1Pt, 'color', [col(pj,:) 1], 'Linewidth', 2);
            pj = 2; p(pj) = plot(ax(a), T.time, T.T2Pt, 'color', [col(pj,:) 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t1 = text_corner(ax(a), ['T1P [volt]'], 1);  
               t1.Color = [col(1,:)];
               t2 = text_corner(ax(a), ['T2P [volt]'], 3);  
               t2.Color = [col(2,:)];
               
            linkaxes(ax, 'x');   

            abc='abcdefghijklmnopqrst';
            for a = 1:(size(ax,1)*size(ax,2))
               text_corner(ax(a), abc(a), 9);
               set(ax(a), 'Xtick', ceil(tl(1)):round(diff(tl)/5):floor(tl(2)));
            end
            


            datetick(ax(a), 'keepticks',  'keeplimits');
            
            unit    = chi_get_unit_name(basedir); % get unit name
            t = text_corner(ax(1), ['unit ' unit], -2);
            
            
      

	%_____________________save pic______________________

	print(gcf,[basedir 'pics' filesep 'temp.png' ],'-dpng','-r200','-painters') 
	savefig(gcf,[basedir 'pics' filesep 'temp.fig' ]) 
end
