function [] = generate_pitot_eps(basedir, do_parallel,  time_limits, low_or_high, save_spec)
%  [] = generate_pitot_eps(basedir, do_parallel,  time_limits, low_or_high, save_spec)
%
%     This function drives the pitot epsilon processing
%
%     
%     INPUT
%        basedir     :  path to instrument folder
%        do_parallel :  parallel processing (deault 1)
%        time_limits :  time_limts for averaging
%        low_or_high :  0 for low frequency (default); 1 for high frequency
%        save_spec   :  Shall I save the spectrogram ? (default 0)
%
%
%   created by: 
%        Johannes Becherer
%        Mon Nov 27 17:18:51 PST 2017

if nargin < 3
   time_limits = [datenum(2000, 1, 1, 0, 0, 0) datenum(2060, 1, 1, 0, 0, 0)];
end
if nargin < 4
   low_or_high = 0;
end
if nargin < 5
   save_spec   = 0;
end


   %_____________________get list of all raw data______________________
      [fids, fdate] = chi_find_rawfiles(basedir);


   %_____________processing loop through all raw files__________________

   disp('Calculating epsilon based on pitot data ')

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               if low_or_high % high frequency estimate
                  proc_pitot_eps(basedir, fids{f}, 2/(24*3600), [2.5 5], save_spec);
               else
                  proc_pitot_eps(basedir, fids{f}, 1/24/12, [.02 .05], save_spec);
               end
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
            if low_or_high % high frequency estimate
               proc_pitot_eps(basedir, fids{f}, 2/(24*3600), [2.5 5], save_spec);
            else
               proc_pitot_eps(basedir, fids{f}, 1/24/12, [.02 .05], save_spec);
            end
         end
      end

   %____________________merge individual files______________________
   if low_or_high % high frequency estimate
      chi_merge_and_avg(basedir, 'pitot_eps2sec', 0, time_limits);
      if save_spec
         chi_merge_and_avg(basedir, 'spec_pitot2sec', 0, time_limits);
      end
   else
      chi_merge_and_avg(basedir, 'pitot_eps300sec', 0, time_limits);;
      if save_spec
         chi_merge_and_avg(basedir, 'spec_pitot300sec', 0, time_limits);
      end
   end
