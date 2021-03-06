function [N2, Sz, s_TS] = cal_N2_from_TS( TSP_time, T, S, P,  Tz_time, Tz,  dt, ChipodDepth)
%%    [N2, Sz, s_TS] = cal_N2_from_TS( TSP_time, T, S, P,  Tz_time, Tz,  dt)
%
%        This function calculates N2 based on a given temperature gradient
%        and a T-S-relation 
%
%        INPUT
%           TSP_time    :  time vector of T,S, and P
%           T           :  temperature 
%           S           :  salinity
%           P           :  pressure in [dbar]
%           Tz_time     :  time vector of the temperature gradient
%           Tz          :  temperature gradient
%           dt          :  timeintrevals of T-S-relation (in sec)
%           ChipodDepth :  depth of chipod (only used if pressure is bad)
%
%        OUTPUT (all output quantities are on the Tz_time grid)
%           N2          :  N2
%           Sz          :  salinity gradient
%           s_TS        :  slope of T-S relation
%
%   created by: 
%        Johannes Becherer
%        Mon Nov 28 18:27:19 PST 2016

%---------------------spli T and S into pieces dt long pieces----------------------
   J{1}  =  1:length(TSP_time);
   Nf    = round( dt/( diff(TSP_time(1:2))*3600*24  ) );   % Nf is the length of the fragment
   I     = split_fragments(J, Nf, 0);  

%_____________________calculated T-S-slope______________________
   time_sl   = nan(1,length(I));
   s_TS_tmp  = nan(1,length(I));
   for i = 1:length(I)
      time_sl(i) = nanmean(TSP_time(I{i}));
      % slope
      if all(isnan(T(I{i})) == 1) ...
              | all(isnan(S(I{i}))) == 1
          p(1) = NaN;
      else
          % use MATLABs centering & scaling transformation to prevent warnings.
          [p, ~, mu] = polyfit( T(I{i}), S(I{i}),1);
          p(1) = p(1)/mu(2); % undo matlab scaling
      end
      s_TS_tmp(i) = p(1);
   end
   % interpolate slope on commen time 
   s_TS  =  interp1( time_sl, s_TS_tmp, Tz_time);

%_____________________cal salinity gradient______________________

Sz =  s_TS .* Tz;
   
%_____________________calculate N2______________________

   % SJW 24-Mar-2020: There's a bug here that if the pressure sensor dies, 
   % P is all NaN when being fed into the lines below that calculate alpha
   % and beta, resulting in those values being NaN along with N2 being NaN. 
   % Fixing this with an if-statement that uses the manually input ChipodDepth.
   if isnan(nanmean(P))
       alpha = nanmean(sw_alpha(S, T,  ChipodDepth, 'temp'));
       beta  = nanmean(sw_beta(S, T,  ChipodDepth, 'temp'));
   else
       alpha = nanmean(sw_alpha(S, T,  P, 'temp'));
       beta  = nanmean(sw_beta(S, T,  P, 'temp'));
   end      
   g     = 9.81;

   N2 = -g*( -alpha*Tz + beta*Sz );

