%     this is script generates header files 
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 14 11:30:20 PST 2017

clear all;
close all;

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   
%_____________________get header files______________________
   % general headers
      hfid = [basedir filesep 'calib' filesep 'header.mat'] ;
      if exist(hfid, 'file')
         disp(['The standard header ' hfid ' exists already!!'])
         load(hfid)
      else % no header found
         choice = questdlg('I did not find a header file', ...
                  'Header file', ...
                   'Use raw-file header (only chipods)','Find header file','Find header file');

         switch choice % What to do?
            case 'Use raw-file header (only chipods)'    
               disp(['generating new header file based on raw_data']);
               [~, head] = raw_load_chipod([basedir filesep 'raw' filesep fids{1}]);
            case 'Find header file'
               [hfile, hdir] = uigetfile('*.*','Load Binary File');
               load([hdir hfile]);
         end
         save(hfid, 'head');
      end

   % Pitot header
      hfid = [basedir filesep 'calib' filesep 'header_p.mat'] ;
      if exist(hfid, 'file')
         disp(['The Pitot header ' hfid ' exists already!!'])
         load(hfid)
      else % no header found
         choice = questdlg('I did not find a Pitot header file', ...
                  'Header file', ...
                   'Type coeffs manually','Find header file','Find header file');

         switch choice % What to do?
            case 'Type coeffs manually'    
               % information box
               uiwait(msgbox({'In order to manually type in the Pitot coeffs' ...
                           'open the png file in the Calibration/Pitot_tubes folder' ...
                           ['use the png file that belongs to Pitot No ' head.sensor_id(7,:)] ...
                            ' ' 'The Calibration folder should be in the deployment directory'}));

               % manual input
               dlg_title = 'Pitot calibration coeffs';
               prompt = {'Slope of Dynamic port:','Slope of Pitch:', 'Slope of Temperature', 'Slope of Pressure'};
               num_lines = 1;
               defaultans = {'0','0', '0', '0'};
               answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

               W.Tilt = [0 str2num(answer{2}) 0 0 0 ];
               W.T    = [0 str2num(answer{3}) 0 0 0 ];
               W.Ps   = [0 str2num(answer{4}) 0 0 0 ];
               % dynamic pressure
               if str2num(answer{1})<1
                  W.Pd = [0 str2num(answer{1}) 0 0 0];
               else
                  W.Pd = [0 1/str2num(answer{1}) 0 0 0];
               end


            case 'Find header file'
               [hfile, hdir] = uigetfile('*.*','Load Binary File');
               load([hdir hfile]);
         end
         save(hfid, 'W');
      end
