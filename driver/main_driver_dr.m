%     this is the driver for the data reduction chi-processing
%
%    created:
%     Johannes Becherer
%     Tue Aug 23 12:26:41 PDT 2016

clear all;
close all;

%_____________________processing flags______________________

   do_parallel =  1;
   do_onboard  =  0;
   do_homework =  1;  


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   



%_____________simulate onboard data reduction__________________

if do_onboard
   % init parallel pool
   if(do_parallel)
      parpool;
      % parallel for-loop
      parfor f=1:length(fids)
            disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
         try % take care if script crashes that the parpoo is shut down
            DR = dr_onboard(basedir, fids{f});

            % save DR
            is1 = strfind(fids{f},'_');  % find under score
            is2 = strfind(fids{f},'.');  % find dot
            savestamp   = [fids{f}((is1):(is2)) 'mat'];
            sdir = [basedir '/proc/DR' filesep];      
            [~, ~, ~] = mkdir(sdir);
            save([sdir 'dr_' savestamp], 'DR');
         catch
            disp(['!!!!!! ' fids{f} ' crashed while processing  !!!!!!' ]);
         end
      end
      % close parpool
      delete(gcp);
   else
      for f=1:length(fids)
         disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
         DR = dr_onboard(basedir, fids{f});

         % save DR
         is1 = strfind(fids{f},'_');  % find under score
         is2 = strfind(fids{f},'.');  % find dot
         savestamp   = [fids{f}((is1):(is2)) 'mat'];
         sdir = [basedir '/proc/DR' filesep];      
         [~, ~, ~] = mkdir(sdir);
         save([sdir 'dr_' savestamp], 'DR');
      end
   end

   %_____________________merge all days______________________
   disp('merge all days')
               ddir = ['DR'];
               chi_merge_and_avg(basedir, ddir, 0);
end

%_____________________post process satelite data______________________
if do_homework
   
   load ../proc/DR.mat;        % get all reduced data
   load ../input/vel_m.mat;    % get velocity array 
      vel = vel_m;
   load ../calib/header.mat;


   disp('doing homework for reduced satelite data....')

   [chi] =  dr_proc(DR, head, vel)


end

