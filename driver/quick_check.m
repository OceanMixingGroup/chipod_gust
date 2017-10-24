%     this is meant to study individual raw files 
%
%   created by: 
%        Johannes Becherer
%        Tue Nov 22 10:21:01 PST 2016

clear all;
close all;

%_____________________flags______________________
ifid             = 1;  % which raw files you would like to look at

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines

%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   
   rfid = fids{ifid};


%_____________________calibrate data______________________
   [data, head, figs] = quick_look([basedir '/raw/'], rfid);

   %print(gcf,['../pics/quick_' rfid '.pdf'],'-dpdf','-painters')
   print(gcf,['../pics/quick_' rfid '.png'],'-dpng','-r100','-painters')
   

%_____________________calibrate pitot including epsilon______________________
   hfid = [basedir filesep 'calib/header_p.mat'];
   if exist(hfid, 'file');
      load(hfid);
      if isfield(W, 'V0')
         [fig, P] = quick_look_pitot(data, W);
          print(gcf,['../pics/quick_pitot_' rfid '.png'],'-dpng','-r100','-painters')
      else
       disp([' V0 does not exist yet. First Run calibrate_pitot']);
      end
   else
      disp([fid ' does not exit']);
   end

