function [P] = make_P_1sec(basedir)
%%    [P] = make_P_1sec(basedir)
%
%  This function calibrates the Pitot tube based on the 1 sec data saved in temp.mat
%
%
%   created by: 
%        Johannes Becherer
%        Thu Mar  8 10:12:05 PST 2018


   [TL]            = whoAmI_timeLimits(basedir);
   time_range      = TL.pitot;

   % if you don't set the basedir function assumes you 
   % in the mfile folder of the instrument directory
   if nargin < 1
      basedir = '../';
   end

   % make sure the path set correctly
   if isempty(which('pitot_calibrate'))
      addpath(genpath([basedir '/mfiles/chipod_gust/software/']));
   end

   f_head = [basedir '/calib/header_p.mat'];
   if exist(f_head)
      load(f_head);
   else
      error('You have to calibrate the Pitot tube first')
   end
   f_temp = [basedir '/proc/temp.mat'];
   if exist(f_temp)
      load(f_temp);
   else
      error('You have have to generate temp.mat first')
   end

   P.time   = T.time;
   P.T   =  zeros(size(P.time));
   if isfield(TL, 'T1')
      if TL.T1< TL.T2
         P.T      = T.T2;
      else
         P.T      = T.T1;
      end
   elseif  isfield(TL, 'T1')
      P.T      = T.T;
   end
   P.depth  =  T.depth;
   P.cmp    = T.cmp;
   P.W      = T.W;
   P.P      = T.P;

   [P.spd, ~, ~] = pitot_calibrate_time( P.time, P.W, P.T, P.P, W);

   ii_nan = (P.spd < 0);

   P.U  = pitot_add_direction( P.time, P.spd, T.time, T.cmp);
      P.spd(ii_nan) =  nan;
      P.U(ii_nan) =  nan;
   P.u  = real(P.U);
   P.v  = imag(P.U);

   save([basedir '/proc/P_1sec.mat'], 'P');
