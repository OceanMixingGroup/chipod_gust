%     this is the driver for Pitot calibration
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 14 11:30:20 PST 2017

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 1;     % use paralelle computing 
   do_raw_data = 1;     % do the averaging of the raw-data (1) or skip (0) if done before
   avg_window = 600;  % averging window for praw (default 600 [sec])
   do_v0_self  = 1;     % detremine V0 based on a min of the averaged signal (self contained)
   DcalWindow =  1000;    % 1000 days window for V0 time (effectivly only a single value for the entire record)
   DcalIncrement = 1000;  % 1000 day increment
   do_v0_adcp  = 1;     % detremin V0 based on a fit against reference velocity (adcp) data
   do_plot     = 0;     % generate some figures in ../pics/ to compare the different velocity estimates
   do_vel_p    = 0;     % which calibration should be used for vel_p (0 none (default), 1: adcp, 2: self)
   do_P1sec    = 0;     % this generates proc/P_1sec.mat based on data from temp.mat

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________set time limits______________________
 % get time limits from whoAmI;
   [TL] =   whoAmI_timeLimits(basedir);
   time_range      = TL.pitot;
 % set manually
   %time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   %time_range(2)  = datenum(2030, 1, 1, 0, 0, 0);

   % calibrate in time range different from valid data time range?
   % if so set limits here just as for time_range.
   % by default, both time ranges are equal.
   cal_time_range = time_range;



   % which temperature sensor to use T1 (1) or if T1 is broken T2 (2) ;  
   % for gusTs (0)
   
   if isfield(TL, 'T1')     % chipod
      % check which sensor longer survives
      if TL.T1(2) < TL.pitot(2) & TL.T2(2)>TL.T1(2)
         use_T = 2;
         % adjust pitot calibration time if temp dies early
         if TL.T2(2) < cal_time_range(2) 
            cal_time_range(2) =  TL.T2(2);
         end
      else
         use_T = 1;
         % adjust pitot calibration time if temp dies early
         if TL.T1(2) < cal_time_range(2) 
            cal_time_range(2) =  TL.T1(2);
         end
      end
   
   else  % gust
     use_T = 0; 
      % adjust pitot calibration time if temp dies early
      if TL.T(2) < cal_time_range(2) 
         cal_time_range(2) =  TL.T(2);
      end
   end

   % manually
   % use_T = 1;


   % shall the pressure calibration be switched off

      if TL.P(2) == TL.master(1)
         use_press = 0;
      else
         use_press = 1;
      end
      % set flag manually 
      %use_press   =  0;
   
      % adjust pitot calibration time if pressure dies early
      if TL.P(2) < cal_time_range(2) & use_press
         cal_time_range(2) =  TL.P(2);
      end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%_____________________do raw data processing______________________
   if do_raw_data
      ischipod = ~isfield(TL,'T');
      generate_praw( basedir, do_parallel, time_range, ischipod, avg_window);
   end
   
%_____________________determine V0______________________
   
   if do_v0_self | do_v0_adcp | do_plot
      determine_v0( basedir, do_v0_self, do_v0_adcp, do_plot, do_vel_p, time_range, use_T, use_press, 'on', DcalWindow, DcalIncrement );
   end

%____________________cal P_1sec.mat_____________________
  if do_P1sec
     make_P_1sec(basedir);
  end

 
