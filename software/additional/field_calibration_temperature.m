function [cal_T]  =  field_calibration_temperature(basedir, rfid, time_ref, temp_ref, Dt, tlim)
%% [cal_T]  =  field_calibration_temperature(basedir, rfid, time_ref, temp_ref, [Dt], [tlim])
%  
%     This function allows for a intercalibration of the temerpature sensor of a gust
%     with filed data. 
%        
%        INPUT
%           basedir  :  basedir of the unit
%           rfid     :  select one of the raw files 
%           teim_ref :  time vector of reference temperature data
%           temp_ref :  reference temperature data
%           tlim     :  time limits to restrict the calibration to (default stat and end of raw_file)
%           Dt       :  time average used for to calibrate the data (default 10min)
%
%        OUTPUT
%           cal_T    :  five element vector with calibration coeffs for head.T
%
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 27 12:58:36 PST 2018


%_____________________read raw file______________________
[raw]   =  raw_load_gust([basedir 'raw/' rfid]);

%---------------------set defaults----------------------
if nargin < 5
   Dt =  1/24/6; % 10 min
end
if nargin < 6
   tlim =  raw.time([1 end]); % raw file limits
end

%_____________________sychronize time steps______________________
dt_raw   =  median(diff(raw.time));
dt_ref   =  median(diff(time_ref));

AW_raw   =  round(Dt/dt_raw);
   if AW_raw < 1
      AW_raw = 1;
   end
AW_ref   =  round(Dt/dt_ref);
   if AW_ref < 1
      AW_ref = 1;
   end

ii_raw   = find( raw.time>= tlim(1) & raw.time<= tlim(2)  );
ii_ref   = find( time_ref>= tlim(1) & time_ref<= tlim(2)  );

%_____________________interpolate raw data on ref_time______________________
raw.T_ref   =  clever_interp( raw.time(ii_raw), movmean(raw.T(ii_raw), AW_raw, 'omitnan'), time_ref(ii_ref) );

p =  polyfit( raw.T_ref, movmean( temp_ref(ii_ref), AW_ref, 'omitnan'), 2);

cal_T =  [p(3) p(2) p(1) 0 0];
