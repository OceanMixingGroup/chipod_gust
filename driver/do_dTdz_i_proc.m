%% This script is meant to do the processing of the  interal dTdz
%   generating ./input/dTdz_i.mat
%
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:19:30 PDT 2017

ticstart = tic;
do_parallel = 1;     % use parallel computing
do_just_merge = 0;   % only merge?
do_proc     = 1;  % do you actually want to do the processing
do_plot     = 1;  % do you want to generate graphical output

dt          = 60;    % sec bits of data for analysis
do_P        = 0;     % use pressure instead of acceleration to get z 
min_dz      = 0.1;   % minimum displacement

% winters d'asaro options
wda_params.do_winters_dasaro = 1; % do sorted gradient required for WDA estimate?
wda_params.wda_dt      = 60; % time chunk over which to average sorted profiles
wda_params.do_P        = do_P; % use pressure sensor instead of accelerometer

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines

%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location

 %_____________________set time limits______________________
 % get time limits from whoAmI;
   [TL] = whoAmI_timeLimits(basedir);
   time_lim = TL.master;

if do_proc
    if ~do_just_merge
        %_____________________get list of all raw data______________________
        [fids, fdate] = chi_find_rawfiles(basedir);

        %_____________processing loop through all raw files__________________
        disp('calculating the internal dTdz');
        % init parallel pool
        if(do_parallel)

            % open the parallel pool
            if do_parallel == 1
               parpool;
            else
               parpool(do_parallel);
            end

            % parallel for-loop
            parfor f=1:length(fids)
                try % take care if script crashes that the parpoo is shut down
                    filedate = datenum(fdate{f}(1:8), 'yymmddhh');
                    if (filedate) >= floor(time_lim(1)) & (filedate) < ceil(time_lim(2))
                        disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
                        chi_generate_dTdz_i(basedir, fids{f}, dt, do_P, min_dz, wda_params);
                    end
                catch ME
                    disp(['!!!!!! ' fids{f} ' crashed while processing  internal dTdz structure !!!!!!' ]);
                    disp(ME);
                end
            end
            % close parpool
            delete(gcp);
        else
            for f=1:length(fids)
                filedate = datenum(fdate{f}(1:8), 'yymmddhh');
                if (filedate) >= floor(time_lim(1)) & (filedate) < (time_lim(2))
                    disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
                    chi_generate_dTdz_i(basedir, fids{f}, dt, do_P, min_dz, wda_params);
                end
            end
        end
    end

    %____________________merge individual files______________________
    chi_merge_and_avg(basedir, 'dTdz', 600);

    if wda_params.do_winters_dasaro
        chi_merge_and_avg(basedir, 'dTdz_w', 0);
    end

    %_____________________cp result to the input directory______________________
    if exist('../proc/dTdz.mat','file') == 2
        ! cp ../proc/dTdz.mat ../input/dTdz_i.mat
    elseif exist('../proc/dTdz_600sec.mat','file') == 2
        ! cp ../proc/dTdz_600sec.mat ../input/dTdz_i.mat
    end

    if exist('../proc/dTdz_w.mat','file') == 2
        ! cp ../proc/dTdz_w.mat ../input/dTdz_w.mat
    end
end

%_____________________plotting______________________
if do_plot
    compare_dTdz
end

disp('Finished processing internal dTdz gradient.')
toc(ticstart)