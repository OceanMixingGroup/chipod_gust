%%  This script runs the chi-processing based on Pitot epsilon in the BBL
%
%  
%   created by: 
%        Johannes Becherer
%        Fri Aug  3 11:53:11 PDT 2018


do_parallel = 1;     % use paralelle computing 
do_raw_proc = 1;     % do the raw data processing 
do_plot     = 1;     % generate a over view plot 


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location

% get time limits from whoAmI;
   [TL] =   whoAmI_timeLimits(basedir);
   time_range      = TL.master;


if do_raw_proc
    save_spec     = 1; % shall the spectrum be saved
    spec_length   =  1; % in seconds 
    frange        =  [1/spec_length 10/spec_length];
    %frange        =  [1 50];

    generate_pitot_eps(basedir, do_parallel,  time_range, spec_length, frange,  save_spec);
end


if do_plot
      
      [fig] = plot_pitot_eps( basedir)
      print(fig,['../pics/Pitot_eps.png'],'-dpng','-r200','-painters')
            
end
