    function [] = proc_eps_chi_proc(basedir, rfid, Eps, S, save_name)
%% [] proc_eps_chi_proc(basedir, rfid, Eps, S, [save_name])
%     This function drives the epsilon processing of the Pitot-tube
%     for a single raw-files 
%
%     INPUT
%        basedir      : unit directory
%        rfid         : raw-file name 
%        S.time   :  time vector of spd 
%        S.spd    :  speed
%        Eps.time :  time vector for epsilon
%        Eps.eps  :  epsilon
%        save_name   :  special label (deafult = '');
%
%   created by: 
%        Johannes Becherer
%        Fri Aug  3 11:34:07 PDT 2018


if nargin < 5
   save_name   =  '';
end

%_____________________preper saving______________________
         is1 = strfind(rfid,'_');  % find under score
         is2 = strfind(rfid,'.');  % find dot
      savestamp   = [rfid((is1):(is2)) 'mat'];
      savedir     = [basedir filesep 'proc' filesep 'eps_chi' save_name filesep];

      % in cases where no 'raw_' preceeds the date for the raw filenames,
      % need different string for savestamp
      if isempty(is1)
         savestamp   = ['_' rfid(1:(is2)) 'mat'];
      end

%_____________________load raw_file______________________
   % load header
   head = chi_get_calibration_coefs(basedir);
   data = chi_calibrate_all([basedir filesep 'raw' filesep rfid], head);

%_____________________for chiipods______________________
   if isfield(data, 'T1')   
       data.T = data.T1;
       data.TP = data.T1Pt;
   end

%_____________________prepare data______________________
   Tp.time = data.time_tp;
   Tp.tp   = data.TPt;
   Tp.spec_floor = data.TP_spec_floor;

   T.time  = data.time;
   T.T     = data.T;
   T.depth = data.depth;

%_____________________process______________________
   [chi, ~]   =  eps_chi_proc(Tp, S, Eps, T);


%---------------------save data----------------------
   [~,~,~] =  mkdir(savedir);
   save([savedir  'eps_chi' save_name savestamp], 'chi');

end
