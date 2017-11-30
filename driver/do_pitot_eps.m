%% this script does all the basic calibration excluding chi-processing and save 1s infomation
%     in ./proc/temp.mat
%  
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:01:26 PDT 2017


do_parallel = 1;     % use paralelle computing 
do_raw_proc = 1;     % do the raw data processing 
do_plot     = 0;     % generate a over view plot 

time_range = [datenum(2000, 1, 1, 0, 0, 0) ...
              datenum(2060, 1, 1, 0, 0, 0)];

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location


if do_raw_proc

    low_or_high = 0; % low frequency estimate
    save_spec   = 0; % shall the spectrum be saved
    generate_pitot_eps(basedir, do_parallel,  time_range, low_or_high, save_spec);

end


if do_plot
      
      spec_length = '300sec';
      [fig] = plot_pitot_eps( basedir, spec_length)
      print(fig,['../pics/Pitot_eps_' spec_length '.png'],'-dpng','-r200','-painters')
            
end
