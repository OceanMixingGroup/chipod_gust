function [] = generate_pitot_eps(basedir, do_parallel,  time_limits, spec_length, frange, save_spec, accfilter)
%  [] = generate_pitot_eps(basedir, do_parallel,  [time_limits, spec_length, frange, save_spec, accfilter])
%
%     This function drives the pitot epsilon processing
%
%     
%     INPUT
%        basedir     :  path to instrument folder
%        do_parallel :  parallel processing (deault 1)
%        time_limits :  time_limts for averaging
%        spec_length :  length of the spectra in seconds
%        frange      :  frequency range to fit in
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
   spec_length   = 2;
end
if nargin < 5
   frange        =  [1/spec_length 10/spec_length];
end
if nargin < 6
   save_spec   = 0;
end
if nargin < 7
   accfilter   = 0;
end


   %_____________________get list of all raw data______________________
      [fids, fdate] = chi_find_rawfiles(basedir);


   %_____________processing loop through all raw files__________________

   disp('Calculating epsilon based on pitot data ')

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
               if accfilter
                  proc_pitot_dissipation(basedir, fids{f}, spec_length/(24*3600), frange, save_spec)
               else
                  proc_pitot_dissipation_noaccfilter(basedir, fids{f}, spec_length/(24*3600), frange, save_spec)
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
            if accfilter
               proc_pitot_dissipation(basedir, fids{f}, spec_length/(24*3600), frange, save_spec)
            else
               proc_pitot_dissipation_noaccfilter(basedir, fids{f}, spec_length/(24*3600), frange, save_spec)
            end
         end
      end

   %____________________merge individual files______________________
      chi_merge_and_avg(basedir, 'pitot_eps', 0, time_limits);
