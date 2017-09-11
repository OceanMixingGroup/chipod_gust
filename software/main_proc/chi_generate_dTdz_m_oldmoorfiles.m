function [Tz_m] = chi_generate_dTdz_m_oldmoorfiles(moor,ChipodDepth,sdir);
% [Tz_m] = chi_generate_dTdz_m_oldmoorfiles(moor,ChipodDepth,sdir);
%
%        This function generates an input file for chi-processing, dTdz_m.m, 
%        saved in directory sdir using mooring files from older deployments.
%        
%        Input
%           moor.time   : time vector
%           moor.depth  : depth in meters from surface at all chipod locations on this mooring
%           moor.u      : zonal velocity (*not used in this code*)
%           moor.v      : meridional velocity (*not used in this code*)
%           moor.dTdz   : temperature gradient at depths of chipods
%           moor.N2     : buoyancy frequency at depths of chipods (May or
%                         may not use salinity. Depends on what variables 
%                         were available when moor structure was created.)
%           ChipodDepth : Depth where THIS chipod was deployed
%           sdir        : directory to save dTdz_m.met to'
%
%           Note: all variables u, v, dTdz, and N2 should have size of 
%               (length(moor.depth),lenth(moor.time))
% 
%        Output
%           Tz_m.time   : time vector  (1min averages)
%           Tz_m.Tz     : temperature gradient
%           Tz_m.N2     : buoyancy freqency
%
%   created by: 
%        Johannes Becherer
%        Fri Sep  2 13:49:20 PDT 2016
%   modified by:
%        Sally Warner to accept mooring files from deployments that were
%        processed with the old processing software



%---------------------create time vector----------------------
   % find beginning
   ts = moor.time(1);
   % find end
   tf = moor.time(end);

   % construct time array
   dt = 1/(24*60);      % one minute
   time = ts:dt:tf;     % time array

   % find delta time for the CTD series
   dtmoor = diff(moor.time(1:2));


%----------extract time series at depth of chipod----------------------
inddepth = find(moor.depth == ChipodDepth);
if isempty(inddepth)
    disp('Error in chi_generate_dTdz_m_oldmoorfiles');
    disp('There needs to be a timeseries at ChipodDepth');
    return
else
    % extract data at the depth of the chipod
    [A,B] = size(moor.dTdz);
    if B>A
        dTdz = moor.dTdz(inddepth,:);
        N2   = moor.N2(inddepth,:);
    else
        dTdz = moor.dTdz(:,inddepth)';
        N2   = moor.N2(:,inddepth)';
    end
end
   

%---------------------cal gradient----------------------
   % for these old deployments, mooring variables were never sampled faster
   % than 1 min, therefore no filtering needed.
       if dtmoor/dt < 1
          dTdz = qbutter(dTdz, dtmoor/dt);
          N2 = qbutter(N2, dtmoor/dt);
       end

   % interpolate time series on common time stamp
   Tz_m.Tz = interp1(moor.time, dTdz, time);
   Tz_m.N2 = interp1(moor.time, N2, time);

   Tz_m.time  = time;

   figure;
   plot(Tz_m.time, Tz_m.N2); hold on;
   plot(Tz_m.time(Tz_m.N2 < 0), Tz_m.N2(Tz_m.N2 < 0));
   plot(xlim, [0, 0], '--', 'color', [1 1 1]*0.6);
   ylabel('N^2')
   title('Negative N^2 in red')
   datetick

%---------------------save----------------------
   save([sdir 'dTdz_m.mat'], 'Tz_m');
