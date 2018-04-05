function [] = chi_generate_dTdz_i(basedir, rfid, varargin)
%% [] = chi_generate_dTdz_i(basedir, rfid, [dt, do_P, min_dz])
%
%     This function calculates vertical stratification
%     based the internal data of the chipod/gusT
%
%     INPUT
%        basedir : unit directory
%        rfid    : raw-file name 
%        dt      : time-window for gradient caluclation in sec (default = 60)
%        do_P    : if 1 pressure is used instead of acceleration for z (default 0)
%        min_dz  : if the standart deviation of z for a given intreval is smaller ...
%                  than min_dz, Tz is set to nan (default min_dz = .1 m)
%
%   created by: 
%        Johannes Becherer
%        Mon Nov 28 16:30:14 PST 2016

%_____________________optional arguments______________________
   if nargin < 3
      dt    = 60; 
      do_P  = 0; 
      min_dz= .1;
   else
      dt = varargin{1};
   end

   if nargin < 4
      do_P  = 0; 
      min_dz= .1;
   else
      do_P = varargin{2};
   end

   if nargin < 5
      min_dz= .1;
   else
      min_dz = varargin{3};
   end

   if nargin < 6
       wda_params.do_winters_dasaro = 1;
       wda_params.wda_dt      = 60; % time chunk over which to average sorted profiles
       wda_params.do_P        = do_P; % use pressure sensor instead of accelerometer
   else
       wda_params = varargin{4};
   end

%_____________________preper saving______________________
         is1 = strfind(rfid,'_');  % find under score
         is2 = strfind(rfid,'.');  % find dot
      savestamp   = [rfid((is1):(is2)) 'mat'];
      savedir     = [basedir filesep 'proc' filesep 'dTdz' filesep];

      % in cases where no 'raw_' preceeds the date for the raw filenames,
      % need different string for savestamp
      if isempty(is1)
         savestamp   = ['_' rfid(1:(is2)) 'mat'];
      end

%_____________________load raw_file______________________
   % load header
   head = chi_get_calibration_coefs(basedir);
   data = chi_calibrate_all([basedir filesep 'raw' filesep rfid], head);

%_____________________time lims______________________
   [TL] = whoAmI_timeLimits(basedir);
   time_lim = TL.master;
   data = apply_time_lims_raw_data(data, TL);

%_____________________pick depth (acc or pressure)______________________
   if do_P % pressure
      data.z = -data.depth;
   else % acceleration
      data.z = -data.a_dis_z;
   end
 
%---------------------split pieces dt long pieces----------------------
   J{1}  =  1:length(data.time);
   Nf    = round( dt/( diff(data.time(1:2))*3600*24  ) );   % Nf is the length of the fragment
   I     = split_fragments(J, Nf, round(Nf*9/10));  

%---------------------run through all segments----------------------
   if isfield(data, 'T') % for gusT
      Tz_i.time   = nan(1,length(I));
      Tz_i.Tz     = nan(1,length(I));
      Tz_i.T      = nan(1,length(I));
      Tz_i.z      = nan(1,length(I));
      for i = 1:length(I)
         Tz_i.time(i) = nanmean(data.time(I{i}));
         Tz_i.T(i)    = nanmean(data.T(I{i}));
         Tz_i.z(i)    = nanmean(data.z(I{i}));

         
         if std(data.z(I{i}))>min_dz ... % check if there is enough vertical variability
                 && ~any(isnan(data.z(I{i}))) && ~any(isnan(data.T(I{i})))
            p = polyfit( data.z(I{i}), data.T(I{i}), 1);
            Tz_i.Tz(i) = p(1);
         else % if there is to lee vertical variation set nan
            Tz_i.Tz(i) = nan;
         end
      end

      Tz_i.S         = ones(length(I),1)*35;
      data.S         = ones(length(data.time),1)*35;
      [Tz_i.N2,Tz_i.Sz,~]  = cal_N2_from_TS(data.time, data.T,  data.S, data.depth, Tz_i.time, Tz_i.Tz, 600);

      %---------------------save data----------------------
      [~,~,~] =  mkdir(savedir);
      save([savedir  'dTdz' savestamp], 'Tz_i');

      % Winters & D'Asaro inference
      % do sorted gradient + save temperature bins
      if wda_params.do_winters_dasaro

          T.time = data.time; T.T = data.T;
          Tp.time = data.time_tp; Tp.tp = data.TPt;

          % average to 1 sec frequency (just like chi estimates)
          ndt = round(1./(diff(data.time(1:2))*86400));
          chi.time = moving_average(T.time, ndt, ndt);
          chi.T = moving_average(T.time, ndt, ndt);

          % get the bins
          Tz_w.wda = do_wda_estimate(wda_params, data, chi, T, Tp);

          % process to get dTdz time series
          wda_proc = process_wda_estimate(chi, Tz_w.wda);

          % save the dTdz time series
          Tz_w.Tz = wda_proc.dTdz;
          Tz_w.zT = wda_proc.dzdT;
          Tz_w.time = wda_proc.time;

          savedir     = [basedir filesep 'proc' filesep 'dTdz_w' filesep];
          [~,~,~] = mkdir(savedir);
          save([savedir  'dTdz_w' savestamp], 'Tz_w');
      end

   else % for chipods
      Tz_i.time   = nan(1,length(I));
      Tz_i.z      = nan(1,length(I));

      Tz_i.Tz1    = nan(1,length(I));
      Tz_i.T1     = nan(1,length(I));
      Tz_i.Tz2    = nan(1,length(I));
      Tz_i.T2     = nan(1,length(I));
      Tz_i.Tz12   = nan(1,length(I));
      Tz_i.T12    = nan(1,length(I));
      for i = 1:length(I)
         Tz_i.time(i) = nanmean(data.time(I{i}));
         Tz_i.z(i)    = nanmean(data.z(I{i}));

         % T1
         Tz_i.T1(i) = nanmean(data.T1(I{i}));
         % T2
         Tz_i.T2(i)  = nanmean(data.T2(I{i}));
         % combo T1 and T2
         Tz_i.T12(i)  = nanmean( .5*(data.T2(I{i}) + data.T1(I{i})) );
         if nanstd(data.z(I{i})) > min_dz && ~any(isnan(data.z(I{i})))
             if ~any(isnan(data.T1(I{i})))
                 p          = polyfit( data.z(I{i}), data.T1(I{i}),1);
                 Tz_i.Tz1(i) = p(1);
             end
             if ~any(isnan(data.T2(I{i})))
                 p          = polyfit( data.z(I{i}), data.T2(I{i}),1);
                 Tz_i.Tz2(i) = p(1);
             end
             if ~any(isnan(data.T1(I{i}))) && ~any(isnan(data.T2(I{i})))
                 p          = polyfit( data.z(I{i}), .5*(data.T2(I{i}) + data.T1(I{i})) , 1);
                 Tz_i.Tz12(i) = p(1);
             end
         else % if there is to lee vertical variation, or nans in data, set nan
            Tz_i.Tz1(i) = nan;
            Tz_i.Tz2(i) = nan;
            Tz_i.Tz12(i) = nan;
         end
      end

      data.S         = ones(length(data.time),1)*35;
      Tz_i.S         = ones(1,length(I))*35;

      % T1
      [Tz_i.N2_1,Tz_i.Sz,~]  = cal_N2_from_TS(data.time, data.T1,  data.S, data.depth, Tz_i.time, Tz_i.Tz1, 600);

      % T2
      [Tz_i.N2_2,~,~]  = cal_N2_from_TS(data.time, data.T2,  data.S, data.depth, Tz_i.time, Tz_i.Tz2, 600);

      % T12
      [Tz_i.N2_12,~,~]  = cal_N2_from_TS(data.time, .5*(data.T2 + data.T1),  data.S, data.depth, Tz_i.time, Tz_i.Tz12, 600);
   

      %---------------------save data----------------------
      [~,~,~] =  mkdir(savedir);
      save([savedir  'dTdz' savestamp], 'Tz_i');

      % Winters & D'Asaro inference
      % do sorted gradient + save temperature bins
      if wda_params.do_winters_dasaro

          T1.time = data.time; T1.T = data.T1;
          Tp1.time = data.time_tp; Tp1.tp = data.T1Pt;
          T2.time = data.time; T2.T = data.T2;
          Tp2.time = data.time_tp; Tp2.tp = data.T2Pt;

          % average to 1 sec frequency (just like chi estimates)
          ndt = round(1./(diff(data.time(1:2))*86400));
          chi1.time = moving_average(T1.time, ndt, ndt);
          chi1.T = moving_average(T1.time, ndt, ndt);
          chi2.time = moving_average(T2.time, ndt, ndt);
          chi2.T = moving_average(T2.time, ndt, ndt);

          % get the bins
          Tz_w.wda1 = do_wda_estimate(wda_params, data, chi1, T1, Tp1);
          Tz_w.wda2 = do_wda_estimate(wda_params, data, chi2, T2, Tp2);

          % process to get dTdz time series
          wda_proc1 = process_wda_estimate(chi1, Tz_w.wda1);
          wda_proc2 = process_wda_estimate(chi2, Tz_w.wda2);

          % get signs
          if exist([basedir filesep 'input' filesep 'dTdz_m.mat'], 'file')
              load([basedir filesep 'input' filesep 'dTdz_m.mat']);

              Tz_w.sgn_moor = get_wda_sign(wda_proc1.time, Tz_m);
          end

          Tz1.time = Tz_i.time; Tz1.Tz = Tz_i.Tz1;
          Tz2.time = Tz_i.time; Tz2.Tz = Tz_i.Tz2;
          Tz_w.sgn_int1 = get_wda_sign(wda_proc1.time, Tz1);
          Tz_w.sgn_int2 = get_wda_sign(wda_proc2.time, Tz2);

          % if we have mooring gradient let's use that sign
          if isfield(Tz_w, 'sgn_moor')
              sgn1 = Tz_w.sgn_moor;
              sgn2 = sgn1;
              Tz_w.sign_used = 'mooring';
          else
              sgn1 = Tz_w.sgn_int1;
              sgn2 = Tz_w.sgn_int2;
              Tz_w.sign_used = 'internal';
          end

          % save the dTdz time series
          Tz_w.Tz1 = wda_proc1.dTdz .* sgn1;
          Tz_w.zT1 = wda_proc1.dzdT .* sgn1;
          Tz_w.Tz2 = wda_proc2.dTdz .* sgn2;
          Tz_w.zT2 = wda_proc2.dzdT .* sgn2;
          Tz_w.time = wda_proc1.time;

          savedir     = [basedir filesep 'proc' filesep 'dTdz_w' filesep];
          [~,~,~] = mkdir(savedir);
          save([savedir  'dTdz_w' savestamp], 'Tz_w');
      end

   end

end
