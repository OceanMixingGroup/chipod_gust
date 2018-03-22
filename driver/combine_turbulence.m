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

%_____________________get all parameters used for combine______________________
   [CP] = default_parameters_combine_turbulence(basedir);
   %_____________________change as much as you like in CP______________________
      % i.e. 
      %CP.time_range        = [ datenum(2000, 1, 1, 0, 0, 0) datenum(2060, 1, 1, 0, 0, 0)];
      %CP.pflag.proc.mmg_ic = 1;
   CP.pflag.master.winters_dasaro = 1;


%_____________________do everything______________________
CP
do_combine_turb(basedir, savedir, CP, do_combine, do_plot, do_mask, save_fig, is_visible)


