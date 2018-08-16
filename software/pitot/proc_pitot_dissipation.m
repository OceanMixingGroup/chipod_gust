    function [] = proc_pitot_dissipation(basedir, rfid, varargin)
%% [] = proc_pitot_disspation(basedir, rfid, [spec_length], [f_range], [save_spec])
%     This function drives the epsilon processing of the Pitot-tube
%     for a single raw-files 
%
%     INPUT
%        basedir      : unit directory
%        rfid         : raw-file name 
%        spec_length  : spectrum length [days]  (default 5min = 1/(24*12)) 
%        save_spec    : shall I save the actual spectrogram (default 0 NO) 
%
%   created by: 
%        Johannes Becherer
%        Tue Jul 31 15:35:32 PDT 2018


%_____________________default parameters______________________
   if nargin < 3
      spec_length = 1/24/3600; 
   else
      spec_length = varargin{1};
   end
   if nargin < 4
      Peps.f_range =  [1/spec_length 10/spec_length];
   else
      Peps.f_range = varargin{2};
   end
   if nargin < 5
      save_spec = 0;
   else
      save_spec = varargin{3};
   end


%_____________________preper saving______________________
         is1 = strfind(rfid,'_');  % find under score
         is2 = strfind(rfid,'.');  % find dot
      savestamp   = [rfid((is1):(is2)) 'mat'];
      savedir     = [basedir filesep 'proc' filesep 'pitot_eps'  filesep];

      % in cases where no 'raw_' preceeds the date for the raw filenames,
      % need different string for savestamp
      if isempty(is1)
         savestamp   = ['_' rfid(1:(is2)) 'mat'];
      end

%_____________________load raw_file______________________
   % load header
   head = chi_get_calibration_coefs(basedir);
   data = chi_calibrate_all([basedir filesep 'raw' filesep rfid], head);

%_____________________calibrate pitot______________________
   fid_pitot_header = [basedir filesep 'calib' filesep 'header_p.mat'];
   if exist(fid_pitot_header, 'file');
      load(fid_pitot_header);
   else
      error([fid_pitot_header  ' does not exit.']);
   end

   if isfield(data, 'T1')   % for chipods
       data.T = data.T1;
   end
   
   [data.spd, ~, ~] = pitot_calibrate_time(data.time, data.W, data.T, data.P, W);


   dt = median( diff(data.time));
   DT = spec_length/2;      
   f_sample = round(1/(dt*3600*24));
      %N_spec = round(DT/dt);
      N_spec = 2^round(log2(round(DT/dt)));
      J{1} = 1:length(data.time);
      [I] = split_fragments( J, N_spec, round(N_spec/2) );

      % find frange
         [~, f] = gappy_psd( data.spd(I{1}) , N_spec, f_sample, 10);
         ii_frange   =  find(f>=Peps.f_range(1) & f<=Peps.f_range(2));


   % initialize fields
   Peps.time = nan(1,length(I));
   Peps.spd  = Peps.time;
   Peps.eps  = Peps.time;
   Peps.eta  = Peps.time;
   Peps.ks_range = nan(2,length(I));
   if save_spec
      %Peps.D_k   = nan(N_spec/2,length(I));
      %Peps.k     = nan(N_spec/2,length(I));
      Peps.D_k   = nan(length(ii_frange),length(I));
      Peps.k     = nan(length(ii_frange),length(I));
   end
   
   nu       = 2.3e-6; % This should be dynamically calculated
   for i = 1:length(Peps.time)
      Peps.time(i)  = mean(data.time(I{i}));
      Peps.spd(i)   = mean(data.spd(I{i}));

      if Peps.spd(i) > 0
         [E, f] = gappy_psd( data.spd(I{i}) , N_spec, f_sample, 10);
            f = f(ii_frange);
            E = E(ii_frange);
         k = 2*pi/Peps.spd(i)*f;
         D = 2*nu*k.^2.*E*Peps.spd(i)/2/pi; 
         [Peps.eps(i), eta, ~] =   dissipation_spec(k, D);
         % if the range is very different repeat calculation
         if k(end)*eta<1e-1  
            [Peps.eps(i), eta, ~] =   dissipation_spec(k, D, eta);
         end
         Peps.eta(i) = eta;
         Peps.ks_range(:,i) = k([1 end])*eta;
         if save_spec
            if ~isfield(Peps, 'f') % f does not change
                Peps.f        = f;
            end
            Peps.D_k(:,i) = D;
            Peps.k(:,i)   = k;
         end
      end
   end


   %_____________________flag data______________________
   Peps = flag_Peps(Peps); 


%---------------------save data----------------------
   [~,~,~] =  mkdir(savedir);
   save([savedir  'pitot_eps' savestamp], 'Peps');

end
