    function [] = proc_pitot_dissipation_noaccfilter(basedir, rfid, varargin)
%% [] = proc_pitot_disspation_noaccfilter(basedir, rfid, [spec_length], [f_range], [save_spec])
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
   DT = spec_length;      
   f_sample = round(1/(dt*3600*24));
   segment_overlap = .5;   % 50 %
      %N_spec = round(DT/dt);
      %N_spec = 2^round(log2(round(DT/dt)));
      N_spec = (round(DT/dt));
      Nfft   = 2^floor(log2(N_spec/3)); 
      J{1} = 1:length(data.time);
      [I] = split_fragments( J, N_spec, round(N_spec*segment_overlap) );

      % find frange
         [~, f] = gappy_psd( data.spd(I{1}) , N_spec, f_sample, 10);
         ii_frange = find(f>=Peps.f_range(1) & f<=Peps.f_range(2));


   % initialize fields
   Peps.time = nan(1,length(I));
   Peps.spd  = Peps.time;
   Peps.T    = Peps.time;
   Peps.depth= Peps.time;
   Peps.eps  = Peps.time;
   Peps.nu   = Peps.time;
   Peps.eta  = Peps.time;
   Peps.ks_range = nan(2,length(I));
   if save_spec
      Peps.D_k_acx   = nan(length(ii_frange),length(I));
      Peps.D_k_acy   = nan(length(ii_frange),length(I));
      Peps.D_k_acz   = nan(length(ii_frange),length(I));
      Peps.D_k   = nan(length(ii_frange),length(I));
      Peps.D_k_org   = nan(length(ii_frange),length(I));
      Peps.k     =  nan(length(ii_frange),length(I));
   end


   
   %nu = 2.3e-6; % This should be dynamically calculated


   for i = 1:length(Peps.time)
      N = length(I{i});
      Peps.time(i)  = mean(data.time(I{i}));
      Peps.spd(i)   = nanmean(data.spd(I{i}));
      Peps.var_acc(i)   = nanvar( movmean(data.a_vel_x(I{i}), round(N/4)))+nanvar( movmean(data.a_vel_y(I{i}), round(N/4)));
      Peps.T(i)     = nanmean(data.T(I{i}));
      Peps.depth(i) = nanmean(data.depth(I{i}));
      Peps.nu(i) = sw_visc(35, Peps.T(i), Peps.depth(i));
   
      f2krad   = 2*pi/Peps.spd(i);
      f2kcyc   = f2krad/(2*pi);

      if Peps.spd(i) > 0.05 
         % pitot 
         tmp_spd     =  data.spd(I{i}) ;
         tmp_spd(tmp_spd<.02) = nan;
         [E11, f]     = gappy_psd( tmp_spd, N_spec, f_sample, 10);
            if length(f)<max(ii_frange)
               continue;
            end
            f = f(ii_frange);
            E11     = E11(ii_frange);

            % transform from u = dudx
            D11f     = E11.*(f2krad*f).^2;

         krad = f2krad*f;
         kcpm = f/Peps.spd(i);
         
         D11k     = D11f/f2krad;
         % I think 2 is wrong ... it should be 7.5 for 1D spec (transferse) and 15 for 
         % longitudinal see pope p134? 
         %D = 2*nu*krad.^2 .*E/f2krad; 

         D2eps = 15*Peps.nu(i);
         %D = D2eps*krad.^2 .*E11/f2krad; 

         [Peps.eps(i), eta, ~] =   dissipation_spec(krad, D11k, Peps.nu(i));
         % if the range is very different repeat calculation
         if krad(end)*eta<1e-1  
            [Peps.eps(i), eta, ~] =   dissipation_spec(krad, D11k, Peps.nu(i), eta);
         end
         Peps.eta(i) = eta;
         Peps.ks_range(:,i) = krad([1 end])*eta;
         if save_spec
            if ~isfield(Peps, 'f') % f does not change
                Peps.f        = f;
            end
            Peps.k(:,i)       = krad;
            Peps.D_k(:,i)     = D11k;
            Peps.D_k_org(:,i) = D11k;
         % acc
            [E11_acx, ~] = gappy_psd( data.a_vel_x(I{i}), N_spec, f_sample, 10);
            [E11_acy, ~] = gappy_psd( data.a_vel_y(I{i}), N_spec, f_sample, 10);
            [E11_acz, ~] = gappy_psd( data.a_vel_z(I{i}), N_spec, f_sample, 10);
            Peps.D_k_acx(:,i) = E11_acx(ii_frange)*f2krad.*(f).^2;
            Peps.D_k_acy(:,i) = E11_acy(ii_frange)*f2krad.*(f).^2;
            Peps.D_k_acz(:,i) = E11_acz(ii_frange)*f2krad.*(f).^2;
         end
      end
   end


   %_____________________flag data______________________

   Peps = flag_Peps(Peps); 


%---------------------save data----------------------
   [~,~,~] =  mkdir(savedir);
   save([savedir  'pitot_eps' savestamp], 'Peps');

end
