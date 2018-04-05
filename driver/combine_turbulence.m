%     this script is meant to combine all processed turbulence information 
%
%    created:
%     Johannes Becherer
%     Tue Aug 23 12:26:41 PDT 2016

clear all;
close all;

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines

%_____________________set flags______________________
   do_combine  =  1; % do actually the averaging calculation (can take a couple o minutes)
   do_plot     =  1; % generate a comparison plot between the different estimates
   do_mask     =  1; % NaN chi estimates using min dTdz, speed thresholds
   save_fig    =  0; % save .fig files?
   is_visible  =  'on'; % shall the figures be displayed

%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data

   filename = [basedir 'mfiles/out/out_combine_turbulence_' datestr(now) '.txt'];
   diary(filename);

%_____________________get all parameters used for combine______________________
   [CP] = default_parameters_combine_turbulence(basedir);

   %========================================================================
   % Check default_parameters_combine_turbulence.m for all available options
   % and comments on how to use them
   %========================================================================

   %_____________________change as much as you like in CP______________________

      %CP.time_range        = [datenum(2000, 1, 1, 0, 0, 0) datenum(2060, 1, 1, 0, 0, 0)];
      %CP.T1death           = datenum(2000, 1, 1, 0, 0, 0);
      %CP.T2death           = datenum(2000, 1, 1, 0, 0, 0);
      %CP.adcpdeath         = datenum(2000, 1, 1, 0, 0, 0);
      %CP.nantimes{1}       = []; % sensor T1 for chipod or T sensor on gusT
      %CP.nantimes{2}       = []; % sensor T2 for chipod
      %CP.nantimes{3}       = []; % pitot sensor
      %CP.avgwindow         = 600;
      %CP.depth             = 0;
      %
      %  switch of master pflags
      %CP.pflag             = CP.pflag.c_ic(0);
      %CP.pflag             = CP.pflag.c_vc(0);
      %CP.pflag             = CP.pflag.c_T1(0);
      %CP.pflag             = CP.pflag.c_T2(0);
      %CP.pflag             = CP.pflag.c_Tzi(0);
      %CP.pflag             = CP.pflag.c_Tzm(0);
      %CP.pflag             = CP.pflag.c_vel_m(0);
      %CP.pflag             = CP.pflag.c_vel_p(0);
      %
      %  or specific flags like
      %CP.pflag.proc.mmg_ic = 1;

   % Tolerance factor for near noise-floor Tp observations
   % CP.factor_spec_floor = 1.5;

   % Fill value for when instrument is at noise floor.
   % think about setting this to 0. NaN recovers previous behaviour
   % CP.noise_floor_fill_value = nan;

   % mask out fits in the IC range with 1 sec data?
   % this tends to happen when dTdz -> 0.
   % Removes a bunch of high values that are suspicious.
   % CP.mask_ic_fits = 1;

   CP.pflag.master.winters_dasaro = 0;


%_____________________do everything______________________
CP
do_combine_turb(basedir, savedir, CP, do_combine, do_plot, do_mask, save_fig, is_visible)
diary off
