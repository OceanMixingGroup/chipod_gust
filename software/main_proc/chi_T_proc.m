function [] = chi_T_proc(basedir, rfid, varargin)
%% [] = chi_T_proc(basedir, rfid, [Dtime])
%     This function drives the processing for the general variables
%     of a single raw-files based on the given flags
%
%     INPUT
%        basedir : unit directory
%        rfid    : raw-file name 
%        Dtime   : optional argument to at a time shift for chipod/gust data
%
%   created by: 
%        Johannes Becherer
%        Wed Sep 21 11:18:51 PDT 2016

%_____________________time shift______________________
   if nargin < 4
      Dtime = 0; 
   else
      Dtime = varargin{1};
   end


%_____________________preper saving______________________
         is1 = strfind(rfid,'_');  % find under score
         is2 = strfind(rfid,'.');  % find dot
      savestamp   = [rfid((is1):(is2)) 'mat'];
      savedir     = [basedir filesep 'proc' filesep 'temp' filesep];

      % in cases where no 'raw_' preceeds the date for the raw filenames,
      % need different string for savestamp
      if isempty(is1)
         savestamp   = ['_' rfid(1:(is2)) 'mat'];
      end

%_____________________load raw_file______________________
   % load header
   head = chi_get_calibration_coefs(basedir);
   data = chi_calibrate_all([basedir filesep 'raw' filesep rfid], head);
%_____________________average on 1 sec intervalls______________________
   T = average_fields(data, 1);

   % cal temperature variance separatly
   dt_p        = median(diff(data.time_tp));
   Np          = round(1/(dt_p*36000*24))*10; % how many valus a second; 
   Np_half     = round(Np/2);
   varT1_tmp   =  movvar( data.T1Pt, Np, 'omitnan');
   T.varT1p    =  varT1_tmp(Np_half:Np:end-Np_half);
   varT2_tmp   =  movvar( data.T2Pt, Np, 'omitnan');
   T.varT2p    =  varT2_tmp(Np_half:Np:end-Np_half);


%---------------------save data----------------------
   [~,~,~] =  mkdir(savedir);
   save([savedir  'temp' savestamp], 'T');


end
