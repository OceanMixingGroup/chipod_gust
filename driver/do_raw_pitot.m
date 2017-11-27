do_parallel = 1;

%_____________________include path of processing files______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines

%____________________set directories______________________
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);

%_____________________average all Pitot raw data to 10 min bins______________________
disp('averaging Pitot raw  data in ./proc/Praw.mat ')
% init parallel pool

generate_praw( basedir, do_parallel)
