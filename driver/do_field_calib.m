%% this script drives field calibration

%_____________________which raw file______________________
i=2;
do_save  =  1; % shall the calibration data be written to the header file?


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir    = [basedir filesep 'raw' filesep]; % raw files location

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);

%_____________________load header ______________________
load([basedir 'calib/header.mat']);

%_____________________load reference temperature______________________
load([basedir 'input/reference_temp.mat']);

%_____________________do field calibration______________________
Dt = 10/24/60; 
[cal_T]  =  field_calibration_temperature(basedir, fids{i}, T.time, T.T, Dt);

if do_save
   disp(['the old coefs are ' num2str(head.coef.T)]);
   disp(['the new coefs are ' num2str(cal_T)]);
   head.coef.T =  cal_T;
   save([basedir 'calib/header.mat'], 'head');
end



