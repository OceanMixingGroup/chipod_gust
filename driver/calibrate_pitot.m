%     this is the driver for Pitot calibration
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 14 11:30:20 PST 2017

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 0;     % use paralelle computing 
   do_raw_data = 0;     % do the averaging of the raw-data (1) or skip (0) if done before
   do_v0_self  = 0;     % detremine V0 based on a min of the averaged signal (self contained)
   do_v0_adcp  = 0;     % detremin V0 based on a fit against reference velocity (adcp) data
   do_plot     = 0;     % generate some figures in ../pics/ to compare the different velocity estimates
   do_vel_p    = 0;     % which calibration should be used for vel_p (0 none (default), 1: adcp, 2: self)

   % This is the time range where the pitot sensor is returning
   % good data
   time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   time_range(2)  = datenum(2030, 1, 1, 0, 0, 0);

   % calibrate in time range different from valid data time range?
   % if so set limits here just as for time_range.
   % by default, both time ranges are equal.
   cal_time_range = time_range;

   % which temperature sensor to use T1 (1) or if T1 is broken T2 (2) ;  
   % for gusTs (0)
   use_T = 1;  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%_____________________include path of processing files______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________do raw data processing______________________
   if do_raw_data
      generate_praw( basedir, do_parallel, time_range);
   end
   
%_____________________determine V0______________________
   
   if do_v0_self | do_v0_adcp | do_plot
      determine_v0( basedir, do_v0_self, do_v0_adcp, do_plot, do_vel_p, time_range, use_T )
   end

 