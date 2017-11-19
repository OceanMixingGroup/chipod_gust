function [] = generate_temp( basedir, do_parallel,  timelims)
%%    [] = generate_temp( basedir,  [do_parallel], [timelims])
%
%     This function generates the temp.mat files in ./proc/
%
%     INPUT
%        basedir     :  path to the directory of the instrument
%        do_parallel :  parallel computing? (default 0 no)
%        timelims    :  optionial argument to restrict the time range
%
%
%   created by: 
%        Johannes Becherer
%        Fri Nov 17 18:37:00 PST 2017

if nargin < 2
   do_parallel = 0;
end

if nargin < 3
   timelims = [datenum(1900,1,1) datenum(2100,0,0)];  % not used yet
end


   %_____________________get list of all raw data______________________
      [fids, fdate] = chi_find_rawfiles(basedir);


   %_____________processing loop through all raw files__________________

   disp('calibrating all base quantities in ./proc/temp.mat ')

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               chi_T_proc(basedir, fids{f});
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
            chi_T_proc(basedir, fids{f});
         end
      end

   %____________________merge individual files______________________
      chi_merge_and_avg(basedir, 'temp', 0, timelims);
      % average 600 sec
      chi_merge_and_avg(basedir, 'temp', 600, timelims);

