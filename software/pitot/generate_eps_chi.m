function [] = generate_eps_chi(basedir, Eps, S, do_parallel,  time_limits, save_name)
%  [] = generate_eps_chi(basedir, Eps, S, do_parallel,  time_limits, [save_name])
%
%     This function drives the pitot epsilon processing
%
%     
%     INPUT
%        basedir     :  path to instrument folder
%        Eps.time    :  time vector for epsilon
%        Eps.eps     :  epsilon
%        S.time      :  time vector of spd 
%        S.spd       :  speed
%        do_parallel :  parallel processing (deault 1)
%        time_limits :  time_limts for averaging
%        save_name   :  special label (deafult = '');
%
%
%   created by: 
%        Johannes Becherer
%        Fri Aug  3 11:45:19 PDT 2018

if nargin < 6
   save_name   =  '';
end
if nargin < 5
   % get time limits from whoAmI;
   [TL]            =   whoAmI_timeLimits(basedir);
   time_range      = TL.master;
end


   %_____________________get list of all raw data______________________
      [fids, fdate] = chi_find_rawfiles(basedir, time_limits);


   %_____________processing loop through all raw files__________________

   disp('Calculating chi based on Tp and epsilon ')

      % init parallel pool
      if(do_parallel)
         if do_parallel > 1
            parpool(do_parallel);
         else
            parpool;
         end
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               proc_eps_chi_proc(basedir, fids{f}, Eps, S, save_name)
            catch ME
               disp(['!!!!!! ' fids{f} ' crashed while processing eps_chi !!!!!!' ]);
               disp(ME)
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            proc_eps_chi_proc(basedir, fids{f}, Eps, S, save_name)
         end
      end

   %____________________merge individual files______________________
      chi_merge_and_avg(basedir, ['eps_chi' save_name], 0, time_limits);
