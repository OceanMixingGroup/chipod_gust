function [CP] = default_parameters_combine_turbulence( basedir)
%% [CP] = default_parameters_combine_turbulence(basedir)
%
%     Sets all parameter necessary for combine turbulence to its default values
%
%     INPUT
%        basedir     :  unit base directory
%     OUT
%        CP          :  structure (combine parameters)
%
%
%   created : 
%        Thu Dec 28 16:59:04 PST 2017
%

   % set thresholds for masking
   CP.min_N2         = 1e-9;
   CP.min_dTdz       = 1e-3;
   CP.min_spd        = 0.05;
   CP.min_inst_spd   = CP.min_spd; % min instantaneous speed past sensor
   CP.mask_inst_spd  = 1; % estimates are crappy if sensor isn't moving
                      % enough.
                      % screws the spectrum calculation...

   % maximum values; anything greater is NaNed out
   CP.max_chi     = 1e-3;
   CP.max_eps     = 1e-2;
   CP.max_Kt      = 1;
   CP.max_Jq      = 1e4;

   CP.avgwindow   = 600; % averaging window in seconds
   CP.avgvalid    = 30; % percent valid values in averaging window for avg to
                  % be non-NaN

   % deglitching parameters
   CP.deglitch_window   = 180; % in seconds
   CP.deglitch_nstd     = 3; % n std. dev. threshold

   % we always mask using speed & dTdz used to calculate chi.
   % the next two are for *additional* masking using a different
   % speed (or dTdz) estimate
   CP.additional_mask_spd  = ''; % '' for none, 'm' for mooring, 'p' for pitot
   CP.additional_mask_dTdz = ''; % '' for none, 'm' for mooring, 'i' for internal

   CP.mask_flushing = 0; % mask so that chipod is always sensing fresh fluid
                      % beta version! turned off by default

   CP.ChipodDepth = 0;

   % normalization for *masking* histograms
   % final processed histograms are always pdf
   CP.normstr = 'count';

   % if you want to restrict the time range that should be combined
   % use the following
   % This restricts time range for BOTH sensors on chipods
   CP.time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   CP.time_range(2)  = datenum(2060, 1, 1, 0, 0, 0);

   % if one sensor dies earlier, specify that time here.
   CP.T1death = datenum(2060, 1, 1, 0, 0, 0); % chipod T1 or gustT T sensor
   CP.T2death = datenum(2060, 1, 1, 0, 0, 0); % T2 sensor



%_________ which estimates should I process?_______________________
   CP.pflag = chi_processing_flags;

      %---------------------gust or chipod----------------------
       CP.pflag = CP.pflag.auto_set(basedir);
      %---------------------add manual flags----------------------
       %CP.pflag = CP.pflag.c_T1(0);       % switch off T1 if bad
       %CP.pflag = CP.pflag.c_T2(0);       % switch off T2 if bad

       %CP.pflag = CP.pflag.c_ic(1);       % switch on ic processing (default off)
       %CP.pflag = CP.pflag.c_vc(0);       % switch off viscous convective processing (default on)
       %CP.pflag.master.epsp = 1;       % switch on eps calculation from pitot (default on)
     
       %CP.pflag = CP.pflag.c_vel_p(0);    % use pitot velocities 
       %CP.pflag = CP.pflag.c_vel_m(0);    % use mooring velocities 
       %CP.pflag = CP.pflag.c_Tzi(0);      % use local (interal) stratification 
       %CP.pflag = CP.pflag.c_Tzm(0);      % use mooring stratification 
      CP.pflag = CP.pflag.make_cons();     % make sub-flags consitent with master flags 



%_____________ additional time ranges to NaN out as necessary________________
   % make an array that looks like
   % nantimes{sensor_number} = [start_time1, end_time1;
   %                            start_time2, end_time2;]
   % start_time & end_time must be datenum
   CP.nantimes{1} = []; % sensor T1 for chipod or T sensor on gusT
   CP.nantimes{2} = []; % sensor T2 for chipod
   CP.nantimes{3} = []; % pitot sensor
