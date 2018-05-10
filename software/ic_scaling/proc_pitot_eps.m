    function [] = proc_pitot_eps(basedir, rfid, varargin)
%% [] = proc_pitot_eps(basedir, rfid, [spec_length], [f_range], [save_spec])
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
%        Wed Sep 21 11:18:51 PDT 2016

%_____________________default parameters______________________
   if nargin < 3
      spec_length = 1/24/12; 
   else
      spec_length = varargin{1};
   end
   if nargin < 4
      Peps.f_range = [.02 .05];
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
      savedir     = [basedir filesep 'proc' filesep 'pitot_eps' num2str(spec_length*24*3600) 'sec' filesep];

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

   [data.spd, Pdym, Vcal] = pitot_calibrate(data.W, data.T, data.P, W);

   %Pdym        = ( data.W - (data.T-W.T0)*W.T(2)  - W.V0 )/W.Pd(2);
   %data.spd    = sign(Pdym).*sqrt(2/1025*abs(Pdym));
   [data.U]    = pitot_add_direction(data.time, data.spd, data.time_cmp, data.cmp);

%___________cal mean speed________________________ 

  dt = median( diff(data.time));
  DT = spec_length/2;      
  data.U_lowpass   = movmean(  data.U , round(DT/dt), 'omitnan' );

   % project in the direction of the mean speed;
   data.U_projected = abs(data.U).*cos(angle(data.U)-angle(data.U_lowpass));
   
   % apply low_pass to wipe out nan smaller than fitting cutoff
    data.U_projected   = movmean(  data.U_projected , round(1/(24*3600*Peps.f_range(2)*4*dt)), 'omitnan' );

%_____________________calculate spectrogram______________________
   [spec, spec_f, Peps.time] = fast_spectrogram(data.time, data.U_projected, spec_length, DT);

%_____________________cal average velocity______________________
   [Peps.vel]   = clever_interp( data.time, data.U, Peps.time );
   [Peps.spd]   = clever_interp( data.time, data.spd, Peps.time );

%_____________________normalize spectrum by ic-scaling______________________
   [ Peps.eps, Peps.var_eps, eps_f] = icscaling_velocity( Peps.time, spec_f, spec, Peps.time, Peps.spd, Peps.f_range);

%_____________________save spectrogram______________________
if save_spec % if you want to save the spectrogram for diagnostics
   Spec.time   =  Peps.time;
   Spec.eps_f  =  eps_f;
   Spec.f      =  spec_f;
   Spec.f_range=  Peps.f_range;
   Spec.eps    =  Peps.eps;
   Spec.U      =  Peps.vel;
   
   spec_dir     = [basedir filesep 'proc' filesep 'spec_pitot'  num2str(spec_length*24*3600) 'sec' filesep];
   [~,~,~] =  mkdir(spec_dir);
   save([spec_dir  'spec_pitot_' num2str(spec_length*24*3600) 'sec' savestamp], 'Spec');;
end

%---------------------save data----------------------
   [~,~,~] =  mkdir(savedir);
   save([savedir  'pitot_eps_' num2str(spec_length*24*3600) 'sec' savestamp], 'Peps');

end
