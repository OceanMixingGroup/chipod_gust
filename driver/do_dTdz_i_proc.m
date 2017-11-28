%% This script is meant to do the processing of the  interal dTdz
%   generating ./input/dTdz_i.mat
%
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:19:30 PDT 2017

do_parallel = 1;     % use paralelle computing 

do_proc     = 1;  % do you actually want to do the processing
do_plot     = 1;  % do you want to generate graphical output

dt          = 60;    % sec bits of data for analysis
do_P        = 0;     % use pressure instead of acceleration to get z 

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location


if do_proc  

 %_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);


 %_____________processing loop through all raw files__________________


   disp('calculating the intrenal dTdz');
   % init parallel pool
   if(do_parallel)
      parpool;
      % parallel for-loop
      parfor f=1:length(fids)
         try % take care if script crashes that the parpoo is shut down
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            chi_generate_dTdz_i(basedir, fids{f}, dt, do_P);
         catch ME
            disp(['!!!!!! ' fids{f} ' crashed while processing  internal dTdz structure !!!!!!' ]);
            disp(ME);
         end
      end
      % close parpool
      delete(gcp);
   else
      for f=1:length(fids)
         disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
         chi_generate_dTdz_i(basedir, fids{f}, dt, do_P);
      end
   end

 %____________________merge individual files______________________
   chi_merge_and_avg(basedir, 'dTdz', 600);

 %_____________________cp result to the input directory______________________
   if exist('../proc/dTdz.mat','file') == 2
      ! cp ../proc/dTdz.mat ../input/dTdz_i.mat
   elseif exist('../proc/dTdz_600sec.mat','file') == 2
      ! cp ../proc/dTdz_600sec.mat ../input/dTdz_i.mat
   end
end

%_____________________plotting______________________
if do_plot
    load ../input/dTdz_i.mat;
       
    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[30 20],'PaperPosition',[0 0 30 20]);
    
            [ax, ~] = create_axes(fig, 3, 1, 0);

            tl = Tz_i.time([1 end]);

            a = 1;
            if isfield(Tz_i,'T')
                plot(ax(a), Tz_i.time, Tz_i.T, 'Linewidth', 1);
                   xlim(ax(a), tl);
            else
                plot(ax(a), Tz_i.time, Tz_i.T1, 'Linewidth', 1);
                   xlim(ax(a), tl);
                plot(ax(a), Tz_i.time, Tz_i.T2, 'Linewidth', 1);
                   xlim(ax(a), tl);  
                plot(ax(a), Tz_i.time, Tz_i.T12, 'Linewidth', 1);
                   xlim(ax(a), tl); 
            end
               t = text_corner(ax(a), ['Temperature [^\circ C]'], 1);
               
            a = 2;
            if isfield(Tz_i,'T')
                po10 = floor(log10(max(abs(Tz_i.Tz))));
                plot(ax(a), Tz_i.time, Tz_i.Tz/10^po10, 'Linewidth', 1);
            else
                po10 = floor(log10(max(abs(Tz_i.Tz1))));
                plot(ax(a), Tz_i.time, Tz_i.Tz1/10^po10, 'Linewidth', 1);
                plot(ax(a), Tz_i.time, Tz_i.Tz2/10^po10, 'Linewidth', 1);
                plot(ax(a), Tz_i.time, Tz_i.Tz12/10^po10, 'Linewidth', 1);                
            end               
                plot(ax(a), tl, [0 0],':k', 'Linewidth', 1);
                xlim(ax(a), tl);
                t = text_corner(ax(a), ['T_z [10^{' num2str(po10) '}K/m]'], 1);
               
            a = 3;
            if isfield(Tz_i,'T')
               po10 = floor(log10(max(abs(Tz_i.N2))));
               plot(ax(a), Tz_i.time, Tz_i.N2/10^po10, 'Linewidth', 1);
            else
               po10 = floor(log10(max(abs(Tz_i.N2_1))));
               plot(ax(a), Tz_i.time, Tz_i.N2_1/10^po10, 'Linewidth', 1);
               plot(ax(a), Tz_i.time, Tz_i.N2_2/10^po10, 'Linewidth', 1);
               plot(ax(a), Tz_i.time, Tz_i.N2_12/10^po10, 'Linewidth', 1);
            end
               plot(ax(a), tl, [0 0],':k', 'Linewidth', 1);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['N^2 [10^{' num2str(po10) '} s^{-2}]'], 1);

            datetick(ax(a), 'keeplimits');
            
            t = text_corner(ax(1), ['T_z of unit ' unit], -2);
            

            print(gcf,'../pics/dTdz_i.png','-dpng','-r200','-painters');
            savefig(fig, '../pics/dTdz_i.fig');
            
            
end
