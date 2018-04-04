%     this is the driver for the main chi-processing
%
%    created:
%     Johannes Becherer
%     Tue Aug 23 12:26:41 PDT 2016

clear all;
close all;

%_____________________check processing flags or process?______________________
   do_dry_run    = 0;
   do_just_merge = 0;


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
   time_lim      = TL.master;

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   

%_____________________set processing flags______________________
   pflag = chi_processing_flags;     % get list of processing flags

   %---------------set processing flags automatically----------
   pflag = pflag.auto_set(basedir);

 % parallel computing
   pflag.master.parallel = 1; % 1 use default number of parallel workers
                               % 0 serial processing
                               % 2...n  use n number of parallel workers


   %---------------------add manual flags----------------------
    %pflag = pflag.c_T1(0);       % switch off T1 if bad
    %pflag = pflag.c_T2(0);       % switch off T2 if bad

    %pflag.master.epsp = 1;       % switch on eps calculation from pitot (default on)
  
    %pflag = pflag.c_vel_p(0);    % use pitot velocities 
    %pflag = pflag.c_vel_m(0);    % use mooring velocities 
    %pflag = pflag.c_Tzi(0);      % use local (interal) stratification 
    %pflag = pflag.c_Tzm(0);      % use mooring stratification 

    %pflag.master.pumped = 0;    % do you use a surface pumped mooring (default = 1)

    pflag.master.use_compass = 1; % if 0, assume chipod vane moves it
                                  % into the flow perfectly

    %pflag.master.use_pres = 0; % set to 0, accelerometers are
                               % working, else differentiate
                               % pressure to get speed past sensor


  % IC-estimate
    %pflag = pflag.c_ic(1);      % switch on ic processing (default off)
    %pflag = pflag.c_vc(0);      % switch off viscous convective processing (default on)
    %pflag.master.ic_dt     = 600; % (in seconds) time window for ic-estimate 
    %pflag.master.ic_frange =  [1/50 1/19]; % fiting range for ic-estimate

    pflag.master.parallel = 1; % 1 use default number of parallel workers
                               % 0 serial processing
                               % 2...n  use n number of parallel workers
    pflag = pflag.make_cons();   % make sub-flags consitent with master flags

   %---------------------get flag status----------------------
   pflag.status();



%_____________main processing________________
if ~do_dry_run   

   do_main_processing( basedir, pflag, time_lim, do_just_merge)

   %_____________________combine all chi_data______________________
   combine_turbulence;

end % ~dry run
