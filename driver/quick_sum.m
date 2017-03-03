%     this routine makes a quick summary of all raw data 
%
%   created by: 
%        Johannes Becherer
%        Mon Oct 24 14:12:06 PDT 2016

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 0;     % use paralelle computing 
   do_qsum     = 0;     % generate qsum.mat 
   do_qsum_pl  = 0;     % generate qsum.mat 


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

   


%%%%%%%%%%%%%%%%%%% genearte quik summary file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_qsum
%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);

   %_____________processing loop through all raw files__________________

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
            try % take care if script crashes that the parpoo is shut down
               chi_qsum_proc(basedir, fids{f});
            catch
               disp(['!!!!!! ' fids{f} ' crashed while processing  !!!!!!' ]);
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
            chi_qsum_proc(basedir, fids{f});
         end
      end

   %_____________________merge individual files______________________
      
      % average 600 sec
      chi_merge_and_avg(basedir, 'qsum', 600);
      chi_merge_and_avg(basedir, 'qsum_r', 600);
end


%---------------------Plot data----------------------
if do_qsum_pl

   chi_qsum_plot(basedir, 'on');
   
end
