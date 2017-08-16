%% This script is meant to do the processing of the  interal dTdz
%   generating ./input/dTdz_i.mat
%
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:19:30 PDT 2017

do_parallel = 1;     % use paralelle computing 

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
! cp ../proc/dTdz.mat ../input/dTdz_i.mat

