%     this is the driver for the main chi-processing
%
%    created:
%     Johannes Becherer
%     Tue Aug 23 12:26:41 PDT 2016

clear all;
close all;

%_____________________check processing flags or process?______________________
   do_dry_run = 0;

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   

%_____________________set processing flags______________________
   pflag = chi_processing_flags;     % get list of processing flags

   %---------------set processing flags automatically----------
   pflag = pflag.auto_set(basedir);

   %---------------------add manual flags----------------------
    %pflag = pflag.c_T1(0);       % switch off T1 if bad
    %pflag = pflag.c_T2(0);       % switch off T2 if bad

    %pflag = pflag.c_ic(1);       % switch on ic processing (default off)
    %pflag = pflag.c_vc(0);      % switch off viscous convective processing (default on)
    %pflag.master.epsp = 1;       % switch on eps calculation from pitot (default on)
  
    %pflag = pflag.c_vel_p(0);    % use pitot velocities 
    %pflag = pflag.c_vel_m(0);    % use mooring velocities 
    %pflag = pflag.c_Tzi(0);      % use local (interal) stratification 
    %pflag = pflag.c_Tzm(0);      % use mooring stratification 

    %pflag.master.pumped = 1;    % do you use a surface pumped mooring (default = 0)

    pflag.master.use_compass = 1; % if 0, assume chipod vane moves it
                                  % into the flow perfectly

    pflag.master.use_pres = 0; % set to 0, accelerometers are
                               % working, else differentiate
                               % pressure to get speed past sensor
    pflag.master.parallel = 1;

    pflag.master.winters_dasaro = 1; % do winters & d'asaro estimate of Kt, Jq
                                     % ONLY for pumped chipods
    pflag.master.wda_dt = 60; % (in seconds) time-interval over which to
                              % apply Winters & D'Asaro methodology

    pflag = pflag.make_cons();   % make sub-flags consitent with master flags

   %---------------------get flag status----------------------
   pflag.status();



%_____________main processing________________
if ~do_dry_run   

   do_main_processing( basedir, pflag)

   %_____________________combine all chi_data______________________
   combine_turbulence;

end % ~dry run
