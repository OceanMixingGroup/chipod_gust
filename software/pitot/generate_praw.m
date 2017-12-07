function [] = generate_praw( basedir, do_parallel,  timelims  , ischipod, avg_window)
%%    [] = generate_praw( basedir,  [do_parallel],  [timelims], [is_chipod], [avg_window])
%
%     This function generates the Praw.mat files in ./proc/
%
%     INPUT
%        basedir     :  path to the directory of the instrument
%        do_parallel :  parallel computing? (default 0 no)
%        timelims    :  optionial argument to restrict the time range
%        is_chipod   :  chipod (1) or gust (0)
%        avg_window  :  averging window for praw (default 600 [sec])
%
%   created by: 
%        Johannes Becherer
%        Mon Nov 27 12:47:43 PST 2017

if nargin < 2
   do_parallel = 0;
end
if nargin < 3
   timelims = [datenum(1900,1,1) datenum(2100,0,0)];  % not used yet
end

if nargin < 4
   % determine if chipod or gusT
   load([ basedir '/calib/header.mat']);
   if isfield(head.coef, 'T1')
       ischipod = 1;
   else
       ischipod = 0;
   end
end

if nargin < 5
   avg_window = 600;  % not used yet
end

   %_____________________get list of all raw data______________________
      [fids, fdate] = chi_find_rawfiles(basedir);


   %_____________processing loop through all raw files__________________
disp('averaging Pitot raw  data in ./proc/Praw.mat ')

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               [~] = pitot_avg_raw_data(basedir, fids{f}, ischipod, avg_window);
            catch ME
               disp(['!!!!!! ' fids{f} ' crashed while processing T structure !!!!!!' ]);
               disp(ME)
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            [~] = pitot_avg_raw_data(basedir, fids{f}, ischipod, avg_window);
         end
      end

   %____________________merge individual files______________________
      chi_merge_and_avg(basedir, 'Praw', 0, timelims);

