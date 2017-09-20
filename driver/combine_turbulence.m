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

   save_fig    = 1; % save .fig files?

   % set thresholds for masking
   min_N2 = 1e-9;
   min_dTdz = 1e-3;
   min_spd = 0.05;
   min_inst_spd = min_spd; % min instantaneous speed past sensor
   mask_inst_spd = 1; % estimates are crappy if sensor isn't moving
                      % enough.
                      % screws the spectrum calculation...

   % maximum values; anything greater is NaNed out
   max_chi = 1e-3;
   max_eps = 1e-2;
   max_Kt = 1;
   max_Jq = 1e4;

   avgwindow = 600; % averaging window in seconds
   avgvalid = 30; % percent valid values in averaging window for avg to
                  % be non-NaN

   % deglitching parameters
   deglitch_window = avgwindow; % in seconds
   deglitch_nstd = 2; % n std. dev. threshold

   % we always mask using speed & dTdz used to calculate chi.
   % the next two are for *additional* masking using a different
   % speed (or dTdz) estimate
   additional_mask_spd = ''; % '' for none, 'm' for mooring, 'p' for pitot
   additional_mask_dTdz = ''; % '' for none, 'm' for mooring, 'i' for internal

   mask_flushing = 0; % mask so that chipod is always sensing fresh fluid
                      % beta version! turned off by default

   ChipodDepth = 30;

   normstr = 'count'; % normalization for histograms

   % if you want to restrict the time range that should be combined
   % use the following
   % This restricts time range for BOTH sensors on chipods
   time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   time_range(2)  = datenum(2030, 1, 1, 0, 0, 0);

   % if one sensor dies earlier, specify that time here.
   T1death = datenum(2060, 1, 1, 0, 0, 0); % chipod T1 or gustT T sensor
   T2death = datenum(2060, 1, 1, 0, 0, 0); % T2 sensor

   % which estimates should I process?
   do = ChooseEstimates(); % get defaults
   % EXAMPLE: turn off all mm estimates
   % do = ChooseEstimates(do, 'no_mm');
   % EXAMPLE: turn off mm1 specifically
   % do.chi_mm1 = 0;

   % additional time ranges to NaN out as necessary
   % make an array that looks like
   % nantimes{sensor_number} = [start_time1, end_time1;
   %                            start_time2, end_time2;]
   % start_time & end_time must be datenum
   nantimes{1} = []; % sensor T1 for chipod or T sensor on gusT
   nantimes{2} = []; % sensor T2 for chipod
   nantimes{3} = []; % pitot sensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

   motionfile = [basedir 'proc' filesep '/motion.mat'];

   % determine if chipod or gusT
   load ../calib/header.mat
   if isfield(head.coef, 'T1')
       isChipod = 1;
   else
       isChipod = 0;
   end

if do_mask
    if additional_mask_dTdz == 'm'
        disp('additional masking using mooring dTdz')
        load([basedir 'input/dTdz_m.mat'])
        Tz = Tz_m;
        clear Tz_m;

    elseif additional_mask_dTdz == 'i'
        disp('additional masking using internal dTdz')
        load([basedir 'input/dTdz_i.mat'])
        Tz = Tz_i;
        clear Tz_i;
    else
        Tz = [];
    end
end

%_____________________find all available chi data______________________
if(do_combine)
   dirname = [savedir 'chi/'];
   d = dir(dirname);
   % eliminate directories
   d = d(~[d(:).isdir]);

   % run through all chi files
   for i = 1:length(d)

      mat_test = strfind(d(i).name, '.mat');
      chi_test = strfind(d(i).name, 'chi');
      ic_test  = strfind(d(i).name, 'ic.');

      if ~isempty(mat_test) & ~isempty(chi_test)

         ID = d(i).name(1:mat_test-1);
         if ~do.(ID)
             disp(['do.' ID ' is disabled!']);
             continue;
         end

         if isChipod
             if ID(end) == '1' | strcmpi(ID(end-3:end), '1_ic') % sensor T1
                 sensor = 1;
             end
             if ID(end) == '2' | strcmpi(ID(end-3:end), '2_ic') % sensor T1
                 sensor = 2;
             end
         else
             sensor = 1; % gusT has only one sensor
         end

         if ~isempty(strfind(ID, 'p'))
             isPitotEstimate = 1;
         else
             isPitotEstimate = 0;
         end

         disp(' ');
         disp(['----------> adding ' ID ]);
         load([dirname ID '.mat'])

         % find desired time range
         iiTrange = find( chi.time >= time_range(1) & chi.time<= time_range(2) );

         % extract valid time range here so that histograms work
         ff = fieldnames(chi);
         for nn = 1:length(ff)
             try
                 chi.(ff{nn}) = chi.(ff{nn})(iiTrange);
             catch ME;
             end
         end

         try
             load ../proc/T_m.mat
             Smean = interp1(T1.time, (T1.S + T2.S)/2, chi.time);
             chi.S = Smean;

             dz = abs(T1.z - T2.z);
             % SBE-37 accuracy is 2e-3 C & 3e-3 psu
             sbe_dTdz = 2*2e-3/dz;
             sbe_dSdz = 2*3e-4/dz;
             sbe_N2 = 9.81 * (1.7e-4 * sbe_dTdz + 7.6e-4 * sbe_dSdz);

             if min_dTdz < sbe_dTdz & ID(6) == 'm'
                 disp(['WARNING: min_dTdz < minimum resolvable based on ' ...
                       'SBE-37 accuracy specifications i.e. ' ...
                       num2str(sbe_dTdz, '%.1e')]);
             end
             if min_N2 < sbe_N2 & ID(6) == 'm'
                 disp(['WARNING: min_N2 < minimum resolvable based on ' ...
                       'SBE-37 accuracy specifications i.e. ' ...
                       num2str(sbe_N2, '%.1e')]);
             end
         catch ME
             Smean = 35*ones(size(chi.time));
         end
         rho = sw_pden(Smean, chi.T, ChipodDepth, 0);
         cp = sw_cp(Smean, chi.T, ChipodDepth);

         %___________________NaN out after sensor death______________________
         if sensor == 1 % sensor T1 or gusT T
             death = find(chi.time > T1death, 1, 'first');
             if ~isempty(death)
                 disp(['NaNing out sensor T1 after it died on ' datestr(T1death)])
                 chi.chi(death:end) = NaN;
                 chi.eps(death:end) = NaN;
             end
         end

         if isChipod
             if sensor == 2 % sensor T2
                 death = find(chi.time > T2death, 1, 'first');
                 if ~isempty(death)
                     disp(['NaNing out sensor T2 after it died on ' datestr(T2death)])
                     chi.chi(death:end) = NaN;
                     chi.eps(death:end) = NaN;
                 end
             end
         end

         %____________________NaN out specific time ranges as necessary____________
         % for temp sensor
         if ~isempty(nantimes{sensor})
             for tt = 1:size(nantimes{sensor}, 1)
                 chi.chi(find_approx(chi.time, nantimes{sensor}(tt, 1), 1): ...
                         find_approx(chi.time, nantimes{sensor}(tt, 2), 1)) = NaN;
             end
             chi.eps(isnan(chi.chi)) = NaN;
         end

         if isPitotEstimate & ~isempty(nantimes{3})
              for tt = 1:size(nantimes{3}, 1)
                 chi.chi(find_approx(chi.time, nantimes{3}(tt, 1), 1): ...
                         find_approx(chi.time, nantimes{3}(tt, 2), 1)) = NaN;
             end
             chi.eps(isnan(chi.chi)) = NaN;
         end

         chi.Kt = 0.5 * chi.chi ./ chi.dTdz.^2;
         chi.Jq = -rho .* cp .* chi.Kt .* chi.dTdz;

         if do_plot
             hfig = CreateFigure;
             Histograms(chi, hfig, normstr, 'raw');

             if ~exist('hfraw', 'var'), hfraw = CreateFigure; end
             Histograms(chi, hfraw, 'pdf', fix_underscore(ID(5:end)));

             if ~exist('hfstrat', 'var')
                 hfstrat = CreateFigure;
                 shown_Tz = '';
             end

             if isempty(strfind(shown_Tz, ID(6)))
                 StratHist(hfstrat, chi, ID);
                 subplot(222);
                 hplt = plot(avgwindow/60*[1, 1], ylim, 'k--');
                 legend(hplt, 'averaging window')
                 shown_Tz = [shown_Tz ID(6)];
             end
         end

         if do_mask
             % speed mask could change depending on estimate
             mask_spd = ID(5);

             if mask_spd == 'm' & ~exist('vel_m', 'var')
                 load ../input/vel_m.mat
                 vel = vel_m;
             elseif mask_spd == 'p' & ~exist('vel_p', 'var')
                 load ../input/vel_p.mat
                 vel = vel_p;
             end
             spdmask = interp1(vel.time, vel.spd, chi.time);

             if additional_mask_spd ~= ''
                 if additional_mask_spd == 'm' & ~exist('vel_m', 'var')
                     load ../input/vel_m.mat
                     vel = vel_m;
                 elseif additional_mask_spd == 'p' & ~exist('vel_p', 'var')
                     load ../input/vel_p.mat
                     vel = vel_p;
                 end
                 addspdmask = interp1(vel.time, vel.spd, chi.time);
             end

             if mask_flushing
                 if ~exist(motionfile, 'file')
                     error(['proc/motion.mat not found. IGNORING ' ...
                            'FLUSHING MASKING. Run do_temp_proc ' ...
                            'to create proc/motion.mat']);
                 end

                 if ~exist('motion', 'var')
                     load(motionfile);
                 end

                 if ~isequal(motion.time, chi.time)
                     ff = fieldnames(motion);
                     for f=1:length(ff)
                         if strcmpi(ff{f}, 'time'), continue; end
                         motion.(ff{f}) = interp1(motion.time, motion.(ff{f}), chi.time);
                     end
                     motion.time = chi.time;
                 end

                 if ~exist('badMotion', 'var') | length(badMotion)~=length(chi.chi)
                     badMotion = make_flushing_mask(motion, mask_spd, vel, do_plot);
                     if do_plot
                         print(gcf, [basedir 'pics' filesep 'angles.png'], '-r200', '-painters', '-bestfit')
                         if save_fig, savefig(gcf, [basedir 'pics' filesep 'angles.fig']); end
                     end
                 end
             end

             if mask_flushing
                 chi = ApplyMask(chi, badMotion, '=', 1, 'volume not being flushed');
                 if do_plot, Histograms(chi, hfig, normstr, 'volume flushed'); end
             end

             [chi, percentage] = ApplyMask(chi, abs(chi.dTdz), '<', min_dTdz, 'Tz');
             chi.stats.dTdz_mask_percentage = percentage;
             perlabel = [' -' num2str(percentage, '%.1f') '%'];
             if do_plot, Histograms(chi, hfig, normstr, ['|Tz| > ' num2str(min_dTdz, '%.1e') perlabel]); end

             [chi, percentage] = ApplyMask(chi, chi.N2, '<', min_N2, 'N2');
             chi.stats.N2_mask_percentage = percentage;
             if percentage > 0.5
                 perlabel = [' -' num2str(percentage, '%.1f') '%'];
                 if do_plot, Histograms(chi, hfig, normstr, ['N2' perlabel]); end
             end

             [chi, percentage] = ApplyMask(chi, chi.spd, '<', min_inst_spd, 'inst speed');
             chi.stats.inst_speed_mask_percentage = percentage;
             [chi, percentage] = ApplyMask(chi, spdmask, '<', min_spd, 'background flow');
             chi.stats.background_flow_mask_percentage = percentage;

             if additional_mask_spd ~= ''
                 [chi, percentage] = ApplyMask(chi, addspdmask, '<', min_spd, 'background flow');
                 chi.stats.additional_background_flow_mask_percentage = percentage;
             end

             % additional Tz masking?
             if ~isempty(Tz)
                 if additional_mask_dTdz == 'i'
                     % choose appropriate internal stratification for sensor
                     Tz.Tz = Tz.(['Tz' ID(7)']);
                 end

                 Tzmask = interp1(Tz.time, Tz.Tz, chi.time);
                 [chi, percentage] = ApplyMask(chi, abs(Tzmask), '<', 1e-4, ...
                                               ['Additional Tz_' additional_mask_dTdz]);
                 chi.stats.additional_dTdz_mask_percentage = percentage;
                 perlabel = [' -' num2str(percentage, '%.1f') '%'];
                 if do_plot, Histograms(chi, hfig, normstr, ['Additional Tz_' additional_mask_dTdz perlabel]); end
             end

             % remove values greater than thresholds
             [chi, chi.stats.max_chi_percentage] = ApplyMask(chi, chi.chi, '>', max_chi, 'max_chi');
             [chi, chi.stats.max_eps_percentage] = ApplyMask(chi, chi.eps, '>', max_eps, 'max_eps');
             [chi, chi.stats.max_Kt_percentage] = ApplyMask(chi, chi.Kt, '>', max_Kt, 'max_Kt');
             [chi, chi.stats.max_Jq_percentage] = ApplyMask(chi, chi.Jq, '>', max_Jq, 'max_Jq');

             if do_plot
                 figure(hfig)
                 set(hfig, 'DefaultLegendBox', 'off');
                 subplot(221); legend(gca, 'show'); title(fix_underscore(ID(5:end)));
                 subplot(222); legend(gca, 'show'); title(fix_underscore(ID(5:end)));
                 subplot(223); legend(gca, 'show'); title(fix_underscore(ID(5:end)));
                 subplot(224); legend(gca, 'show'); title(fix_underscore(ID(5:end)));

                 print(gcf,['../pics/histograms-masking-' ID '.png'],'-dpng','-r200','-painters')
                 if save_fig, savefig(gcf,['../pics/histograms-masking-' ID '.fig']); end
             end
         end

         % convert averaging window from seconds to points
         ww =  round(avgwindow/(diff(chi.time(1:2))*3600*24));
         dw = round(deglitch_window/(diff(chi.time(1:2))*3600*24));

         if isempty(ic_test)
             % deglitch chi and eps before
             % calculating Jq and Kt
             % not required for IC estimate because that is already
             % an averaged estimate
             disp('Deglitch... itch... tch... ch')
             tic;
             chi.chi = deglitch(chi.chi, dw, deglitch_nstd, 'b');
             chi.eps = deglitch(chi.eps, dw, deglitch_nstd, 'b');
             toc;

             % get list of all fields to average
             ff = fields(chi);

             %% average data
             disp('Running moving average')
             tic;
             for f = 1:length(ff)  % run through all fields in chi
                 if ( length(chi.(ff{f})) == length(chi.time) )
                     if strcmpi(ff{f}, 'Kt') | strcmpi(ff{f}, 'Jq'), continue; end
                     Turb.(ID).(ff{f}) = moving_average( chi.(ff{f}), ww, ww , avgvalid);
                 else
                     Turb.(ID).(ff{f}) = chi.(ff{f});
                 end
             end
             toc;

             % reapply dTdz filter
             % (because we divide by this dTdz for Kt, Jq)
             Turb.(ID) = ApplyMask(Turb.(ID), abs(Turb.(ID).dTdz), '<', min_dTdz, 'avg dTdz')

             % recalculate using averaged quantities
             % if we average over a time period greater than
             % sampling period of dTdz, this estimate will differ!
             Turb.(ID).Kt = 0.5 * Turb.(ID).chi ./ Turb.(ID).dTdz.^2 + ...
                 sw_tdif(interp1(chi.time, Smean, Turb.(ID).time), ...
                         Turb.(ID).T, ChipodDepth);
             Turb.(ID).Jq = -1025 .* 4200 .* Turb.(ID).Kt .* Turb.(ID).dTdz;
         else
             Turb.(ID) = chi;
         end

         if do_plot
             if ~exist('hfig2', 'var'), hfig2 = CreateFigure; end
             Histograms(Turb.(ID), hfig2, 'count', fix_underscore(ID(5:end)));
         end
         
         % include statistics
         Turb.(ID).stats.chimean = nanmean(Turb.(ID).chi);
         Turb.(ID).stats.epsmean = nanmean(Turb.(ID).eps);
         Turb.(ID).stats.Ktmean  = nanmean(Turb.(ID).Kt);
         Turb.(ID).stats.Jqmean  = nanmean(Turb.(ID).Jq);
         Turb.(ID).stats.chimedian = nanmedian(Turb.(ID).chi);
         Turb.(ID).stats.epsmedian = nanmedian(Turb.(ID).eps);
         Turb.(ID).stats.Ktmedian  = nanmedian(Turb.(ID).Kt);
         Turb.(ID).stats.Jqmedian  = nanmedian(Turb.(ID).Jq);
         Turb.(ID).stats.chistd = nanstd(Turb.(ID).chi);
         Turb.(ID).stats.epsstd = nanstd(Turb.(ID).eps);
         Turb.(ID).stats.Ktstd  = nanstd(Turb.(ID).Kt);
         Turb.(ID).stats.Jqstd  = nanstd(Turb.(ID).Jq);
      end
   end

   if do_plot
       figure(hfig2)
       subplot(221); title(['Final ' num2str(avgwindow/60) ' min mean']);
       subplot(222); title(['Final ' num2str(avgwindow/60) ' min mean']);

       print(gcf,['../pics/histograms-final.png'],'-dpng','-r200','-painters')
       if save_fig, savefig(gcf,['../pics/histograms-final.fig']); end

       figure(hfraw)
       subplot(221); title(['raw 1s estimates']);
       subplot(222); title(['raw 1s estimates']);
       print(gcf,['../pics/histograms-raw.png'],'-dpng','-r200','-painters')

       figure(hfstrat)
       print(gcf,['../pics/histograms-stratification.png'],'-dpng','-r200','-painters')
   end

   Turb.hash = githash('driver/combine_turbulence.m');
   Turb.do_mask = do_mask;
   Turb.additional_mask_dTdz = additional_mask_dTdz;
   Turb.additional_mask_spd = additional_mask_spd;
   Turb.min_dTdz = min_dTdz;
   Turb.min_spd = min_spd;
   Turb.avgwindow = avgwindow;
   Turb.min_inst_spd = min_inst_spd;
   Turb.min_N2 = min_N2;

   %---------------------add readme----------------------
   Turb.readme = {...
         '------------------sub-fields--------------------'; ...
         'time    :  matlab time'; ...
         'chi     :  chi (10 min averaged deglitched)'; ...
         'eps     :  dissipation (10 min averaged deglitched)'; ...
         'dTdz    :  vertical temperature gradient used for processing'; ...
         'N2      :  buoyancy frequency used for processing'; ...
         'spd     :  flow pass sensor (for pumped moorings contains ACC-information)'; ...
         'T       :  temperature used for calculation'; ...
         'S       :  salinity used for calculation (default 35)'; ...
         'mask    :  values that got flagged during processing due to minTz or speed'; ...
         ''; ...
         '------------------main-fields--------------------'; ...
         '..._ic  :  inertial convective range estimate (no subscript viscous convective range)'; ...
         'mm1     :  mooring speed  - mooring N2              - temp_sensor 1 '; ...
         'mm2     :  mooring speed  - mooring N2              - temp_sensor 2 '; ...
         'mi11    :  mooring speed  - internal N2 (T1 sensor) - temp_sensor 1 '; ...
         'mi112   :  mooring speed  - internal N2 (combination T1 and T2 sensor) - temp_sensor 1 '; ...
         'mi22    :  mooring speed  - internal N2 (T2 sensor) - temp_sensor 2 '; ...
         'pm1     :  Pitot speed    - mooring N2              - temp_sensor 1 '; ...
         'pm2     :  Pitot speed    - mooring N2              - temp_sensor 2 '; ...
         'pi11    :  Pitot speed    - internal N2 (T1 sensor) - temp_sensor 1 ' ...
         };

%_____________________save combined structure______________________
   save([savedir '/Turb.mat'], 'Turb');
end

%_____________________comparison plot______________________
if do_plot    
    % moved this plotting function to a separate function so it can be run
    % independently if needed.
    plot_Turb  
end


% Examples of using TestMask and DebugPlots to check masking
% TestMask(chi, abs(chi.dTdz), '<', [1e-4, 3e-4, 1e-3], 'Tz');
% t0 = datenum(2014, 01, 01);
% t1 = datenum(2014, 03, 01);
% DebugPlots([], t0, t1, chi, 'raw', 1)
% chi1 = ApplyMask(chi, abs(chi.dTdz), '<', 3e-4, 'T_z');
% DebugPlots([], t0, t1, chi1, 'T_z > 3e-4', 1)
% load ../proc/temp.mat
% load ../proc/Turb.mat
