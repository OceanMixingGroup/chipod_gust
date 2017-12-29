function [] = do_combine_turb( basedir, savedir, CP, do_combine, do_plot, do_mask, save_fig, is_visible) 
%%    [] = do_combine_turb( basedir, savedir, CP, do_combine, do_plot, do_mask, save_fig, is_visible) 
%
%        All the fuction of combine turbulence in one routine
%
%        INPUT
%           basedir     :  unit directory
%           savedir     :  where shall the structure be saved
%           CP          :  flagging parameters generated with default_parameters_combine_turbulence.m
%           do_combine  :  1: do_combining 0:don't
%           do_plot     :  generates plots (1 yes 0 no)
%           do_mask     :  apply mask
%           save_fig    :  save fig files
%           is_visible  :  visibilty of figures 'on' or 'off'
%
%
%   created: 
%        Thu Dec 28 18:18:54 PST 2017

%_____________________defaults______________________
if nargin < 8
   is_visible = 'on';
end
if nargin < 7
   save_fig   = 0;
end
if nargin < 6
   do_mask = 1;
end
if nargin < 5
   do_plot  =  1;
end
if nargin < 4
   do_combine;
end
if nargin < 3
   [CP] = default_parameters_combine_turbulence(basedir);
end
if nargin < 2
   savedir = [basdir 'proc/'];
end


   motionfile = [basedir 'proc' filesep '/motion.mat'];

   % determine if chipod or gusT
   load([basedir '/calib/header.mat'])
   if isfield(head.coef, 'T1')
       isChipod = 1;
   else
       isChipod = 0;
   end

if do_mask
    if CP.additional_mask_dTdz == 'm'
        disp('additional masking using mooring dTdz')
        load([basedir 'input/dTdz_m.mat'])
        Tz = Tz_m;
        clear Tz_m;

    elseif CP.additional_mask_dTdz == 'i'
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

         ID = d(i).name(5:mat_test-1);
         if ~CP.pflag.proc.(ID)
             disp([ ID ' is disabled!']);
             continue;
         end

         if isChipod
             if ID(end) == '1' ...
                         | (length(ID) > 3 & strcmpi(ID(end-3:end), '1_ic')) % sensor T1
                 sensor = 1;
             end
             if ID(end) == '2' ...
                         | (length(ID) > 3 & strcmpi(ID(end-3:end), '2_ic')) % sensor T2
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
         load([dirname 'chi_' ID '.mat'])

         % find desired time range
         iiTrange = find( chi.time >= CP.time_range(1) & chi.time<= CP.time_range(2) );

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

             if CP.min_dTdz < sbe_dTdz & ID(2) == 'm'
                 disp(['WARNING: min_dTdz < minimum resolvable based on ' ...
                       'SBE-37 accuracy specifications i.e. ' ...
                       num2str(sbe_dTdz, '%.1e')]);
             end
             if CP.min_N2 < sbe_N2 & ID(2) == 'm'
                 disp(['WARNING: min_N2 < minimum resolvable based on ' ...
                       'SBE-37 accuracy specifications i.e. ' ...
                       num2str(sbe_N2, '%.1e')]);
             end
         catch ME
             Smean = 35*ones(size(chi.time));
         end

         %___________________NaN out after sensor death______________________
         if sensor == 1 % sensor T1 or gusT T
             death = find(chi.time > CP.T1death, 1, 'first');
             if ~isempty(death)
                 disp(['NaNing out sensor T1 after it died on ' datestr(CP.T1death)])
                 chi.chi(death:end) = NaN;
                 chi.eps(death:end) = NaN;
                 chi.T(death:end) = NaN;
             end
         end

         if isChipod
             if sensor == 2 % sensor T2
                 death = find(chi.time > CP.T2death, 1, 'first');
                 if ~isempty(death)
                     disp(['NaNing out sensor T2 after it died on ' datestr(CP.T2death)])
                     chi.chi(death:end) = NaN;
                     chi.eps(death:end) = NaN;
                     chi.T(death:end) = NaN;
                 end
             end
         end

         %____________________NaN out specific time ranges as necessary____________
         % for temp sensor
         if ~isempty(CP.nantimes{sensor})
             for tt = 1:size(CP.nantimes{sensor}, 1)
                 chi.chi(find_approx(chi.time, CP.nantimes{sensor}(tt, 1), 1): ...
                         find_approx(chi.time, CP.nantimes{sensor}(tt, 2), 1)) = NaN;
             end
             chi.eps(isnan(chi.chi)) = NaN;
         end

         if isPitotEstimate & ~isempty(CP.nantimes{3})
              for tt = 1:size(CP.nantimes{3}, 1)
                 chi.chi(find_approx(chi.time, CP.nantimes{3}(tt, 1), 1): ...
                         find_approx(chi.time, CP.nantimes{3}(tt, 2), 1)) = NaN;
             end
             chi.eps(isnan(chi.chi)) = NaN;
         end

         chi.Kt = 0.5 * chi.chi ./ chi.dTdz.^2;
         chi.Jq = -1025 .* 4200 .* chi.Kt .* chi.dTdz;

         % save unmasked chi structure for later use
         chiold = chi;

         if do_plot
             hfig = CreateFigure(is_visible);
             Histograms(chi, hfig, CP.normstr, ID, 'raw');

             if ~exist('hfraw', 'var'), hfraw = CreateFigure(is_visible); end
             Histograms(chi, hfraw, 'pdf', ID, ID);

             if ~exist('hfstrat', 'var')
                 hfstrat = CreateFigure(is_visible);
                 shown_Tz = '';
             end

             if isempty(strfind(shown_Tz, ID(2)))
                 StratHist(hfstrat, chi, ID);
                 subplot(222);
                 hplt = plot(CP.avgwindow/60*[1, 1], ylim, 'k--');
                 legend(hplt, 'averaging window')
                 shown_Tz = [shown_Tz ID(2)];
             end
         end

         if do_mask
             % speed mask could change depending on estimate
             mask_spd = ID(1);

             if mask_spd == 'm' & ~exist('vel_m', 'var')
                 load ../input/vel_m.mat
                 vel = vel_m;
             elseif mask_spd == 'p' & ~exist('vel_p', 'var')
                 load ../input/vel_p.mat
                 vel = vel_p;
             end
             spdmask = interp1(vel.time, vel.spd, chi.time);

             if CP.additional_mask_spd ~= ''
                 if CP.additional_mask_spd == 'm' & ~exist('vel_m', 'var')
                     load ../input/vel_m.mat
                     vel = vel_m;
                 elseif CP.additional_mask_spd == 'p' & ~exist('vel_p', 'var')
                     load ../input/vel_p.mat
                     vel = vel_p;
                 end
                 addspdmask = interp1(vel.time, vel.spd, chi.time);
             end

             if CP.mask_flushing
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
                         if save_fig, set(gcf, 'visible', 'on'); savefig(gcf, [basedir 'pics' filesep 'angles.fig']); end
                     end
                 end
             end

             if CP.mask_flushing
                 chi = ApplyMask(chi, badMotion, '=', 1, 'volume not being flushed');
                 if do_plot, Histograms(chi, hfig, CP.normstr, (ID), 'volume flushed'); end
             end

             [chi, percentage] = ApplyMask(chi, abs(chi.dTdz), '<', CP.min_dTdz, 'Tz');
             chi.stats.dTdz_mask_percentage = percentage;
             perlabel = [' -' num2str(percentage, '%.1f') '%'];
             if do_plot, Histograms(chi, hfig, CP.normstr, (ID), ['|Tz| > ' num2str(CP.min_dTdz, '%.1e') perlabel]); end

             [chi, percentage] = ApplyMask(chi, chi.N2, '<', CP.min_N2, 'N2');
             chi.stats.N2_mask_percentage = percentage;
             if percentage > 0.5
                 perlabel = [' -' num2str(percentage, '%.1f') '%'];
                 if do_plot, Histograms(chi, hfig, CP.normstr,(ID), ['N2' perlabel]); end
             end

             [chi, percentage] = ApplyMask(chi, chi.spd, '<', CP.min_inst_spd, 'inst speed');
             chi.stats.inst_speed_mask_percentage = percentage;
             [chi, percentage] = ApplyMask(chi, spdmask, '<', CP.min_spd, 'background flow');
             chi.stats.background_flow_mask_percentage = percentage;

             if CP.additional_mask_spd ~= ''
                 [chi, percentage] = ApplyMask(chi, addspdmask, '<', CP.min_spd, 'background flow');
                 chi.stats.additional_background_flow_mask_percentage = percentage;
             end

             % additional Tz masking?
             if ~isempty(Tz)
                 if CP.additional_mask_dTdz == 'i'
                     % choose appropriate internal stratification for sensor
                     Tz.Tz = Tz.(['Tz' ID(3)']);
                 end

                 Tzmask = interp1(Tz.time, Tz.Tz, chi.time);
                 [chi, percentage] = ApplyMask(chi, abs(Tzmask), '<', 1e-4, ...
                                               ['Additional Tz_' CP.additional_mask_dTdz]);
                 chi.stats.additional_dTdz_mask_percentage = percentage;
                 perlabel = [' -' num2str(percentage, '%.1f') '%'];
                 if do_plot, Histograms(chi, hfig, CP.normstr, (ID), ['Additional Tz_' CP.additional_mask_dTdz perlabel]); end
             end

             % remove values greater than thresholds
             [chi, chi.stats.max_chi_percentage] = ApplyMask(chi, chi.chi, '>', CP.max_chi, 'max_chi');
             [chi, chi.stats.max_eps_percentage] = ApplyMask(chi, chi.eps, '>', CP.max_eps, 'max_eps');

             if do_plot
                 %figure(hfig)
                 set(0, 'currentfigure', hfig);
                 set(hfig, 'DefaultLegendBox', 'off');
                 subplot(221); legend(gca, 'show'); title((ID));
                 subplot(222); legend(gca, 'show'); title((ID));
                 subplot(223); legend(gca, 'show'); title((ID));
                 subplot(224); legend(gca, 'show'); title((ID));

                 print(gcf,['../pics/histograms-masking-' ID '.png'],'-dpng','-r200','-painters')
                 if save_fig,set(gcf, 'visible', 'on'); savefig(gcf,['../pics/histograms-masking-' ID '.fig']); end
             end
         end

         % convert averaging window from seconds to points
         ww =  round(CP.avgwindow/(diff(chi.time(1:2))*3600*24));
         dw = round(CP.deglitch_window/(diff(chi.time(1:2))*3600*24));

         if isempty(ic_test)
             % deglitch chi and eps before
             % calculating Jq and Kt
             % not required for IC estimate because that is already
             % an averaged estimate
             disp('Deglitch... itch... tch... ch')
             tic;
             chi.chi = 10.^deglitch(log10(chi.chi), dw, CP.deglitch_nstd, 'b');
             chi.eps = 10.^deglitch(log10(chi.eps), dw, CP.deglitch_nstd, 'b');
             toc;

             % get list of all fields to average
             ff = fields(chi);

             %% average data
             disp('Running moving average')
             tic;
             for f = 1:length(ff)  % run through all fields in chi
                 if ( length(chi.(ff{f})) == length(chi.time) )
                     if strcmpi(ff{f}, 'Kt') | strcmpi(ff{f}, 'Jq'), continue; end
                     Turb.(ID).(ff{f}) = moving_average( chi.(ff{f}), ww, ww , CP.avgvalid);
                 else
                     Turb.(ID).(ff{f}) = chi.(ff{f});
                 end
             end
             toc;
         else
             Turb.(ID) = chi;
         end

         % reapply dTdz filter
         % (because we divide by this dTdz for Kt, Jq)
         Turb.(ID) = ApplyMask(Turb.(ID), abs(Turb.(ID).dTdz), '<', CP.min_dTdz, 'avg dTdz');

         % recalculate using averaged quantities
         % if we average over a time period greater than
         % sampling period of dTdz, this estimate will differ!
         Turb.(ID).Kt = 0.5 * Turb.(ID).chi ./ Turb.(ID).dTdz.^2 + ...
             sw_tdif(interp1(chi.time, Smean, Turb.(ID).time), ...
                     Turb.(ID).T, CP.ChipodDepth);
         Turb.(ID).Jq = -1025 .* 4200 .* Turb.(ID).Kt .* Turb.(ID).dTdz;

         [Turb.(ID), Turb.(ID).stats.max_Kt_percentage] = ApplyMask(Turb.(ID), Turb.(ID).Kt, '>', CP.max_Kt, 'max_Kt');
         [Turb.(ID), Turb.(ID).stats.max_Jq_percentage] = ApplyMask(Turb.(ID), Turb.(ID).Jq, '>', CP.max_Jq, 'max_Jq');

         if do_plot
             if ~exist('hfig2', 'var'), hfig2 = CreateFigure(is_visible); end
             Histograms(Turb.(ID), hfig2, 'pdf', ID, ID);
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
       %figure(hfig2)
       set(0, 'currentfigure', hfig2);
       subplot(221); title(['Final ' num2str(CP.avgwindow/60) ' min mean']);
       subplot(222); title(['Final ' num2str(CP.avgwindow/60) ' min mean']);

       print(gcf,['../pics/histograms-final.png'],'-dpng','-r200','-painters')
       if save_fig,set(gcf, 'visible', 'on'); savefig(gcf,['../pics/histograms-final.fig']); end

       %figure(hfraw)
       set(0, 'currentfigure', hfraw);
       subplot(221); title(['raw 1s estimates']);
       subplot(222); title(['raw 1s estimates']);
       print(gcf,['../pics/histograms-raw.png'],'-dpng','-r200','-painters')

       %figure(hfstrat)
       set(0, 'currentfigure', hfstrat);
       print(gcf,['../pics/histograms-stratification.png'],'-dpng','-r200','-painters')
   end

   Turb.hash = githash('driver/combine_turbulence.m');
   Turb.do_mask = do_mask;
   Turb.parameters = CP;

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
    [fig] = plot_Turb(basedir, CP.pflag, is_visible);
      print(fig,[basedir '/pics/Compare_Turb.png'],'-dpng','-r200')
      if save_fig
         disp('saving Compare_Turb.fig')
         set(gcf, 'visible', 'on');
         savefig(fig,[basedir '/pics/Compare_Turb.fig'])
      end
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
end
