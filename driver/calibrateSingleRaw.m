clear all;
close all;

addpath(genpath('./chipod_gust/software/'));


here    =   pwd;                % mfiles folder
basedir =   here(1:(end-6));    % substract the mfile folder

[fids, fdate] = chi_find_rawfiles(basedir);
i=5;

[data, head] = quick_calibrate([basedir 'raw/'], fids{i});

T.time       = data.time;
if isfield(data, 'T')
   T.T          = data.T;
   T.Tp         = data.TPt;
else
   T.T1         = data.T1;
   T.T2         = data.T2;
   T.Tp1        = data.T1Pt;
   T.Tp2        = data.T2Pt;
end

T.P          = data.P;
T.W          = data.W;
T.AX         = data.AX;
T.AY         = data.AY;
T.AZ         = data.AZ;
T.cmp        = data.cmp;
T.time_cmp   = data.time_cmp;



savefid =   [basedir 'raw/mat/cal_' fdate{i} '.mat'];  % directory directory to save data
save(savefid, 'T');
