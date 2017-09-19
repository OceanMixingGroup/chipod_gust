function [Praw] = pitot_avg_raw_data(basedir, rfid, varargin )
%%    [Praw] =  pitot_avg_raw_data(basedir, rfid, [is_chipod, [dt] ])
%
%        This function averages the Pitot tube relevant raw data
%
%  INPUT
%     basedir  :  base directory
%     rfid     :  raw file name
%     is_chipod:  1 for chipod; 0 for gusTs (default 1)
%     dt       :  averaging interval (default 600 sec)
%
%  OUTPUT
%     Praw.time   :  time vector
%     Praw.T1     :  mean raw T1 sensor
%     Praw.T2     :  mean raw T2 sensor
%     Praw.vT1    :  variance of T1 sensor
%     Praw.vT2    :  variance of T2 sensor
%     Praw.P      :  mean pressure
%     Praw.W      :  mean Pitot voltage
%     Praw.cmp    :  mean compass direction
%
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 14 15:17:40 PST 2017





%_____________________parameters______________________

% chipod or gusT default chipod
if nargin < 3                  
   is_chipod = 1;
else
   is_chipod = varargin{1};
end

% what is the averageing window 600 sec default
if nargin < 4                  
   dt = 600;
else
   dt = varargin{2};
end


%_____________________read raw_data______________________
if is_chipod % chipod
   [rdat, ~] = raw_load_chipod([basedir '/raw/', rfid]);
else % gusT
   [rdat] = raw_load_gust([basedir '/raw/', rfid]);
end


%_____________________split data vector______________________
   % time vector
      if is_chipod
         di_cmp        = round(length(rdat.datenum)/length(rdat.CMP)); 
         di            = round(length(rdat.datenum)/length(rdat.T1)); 
         rdat.time     = rdat.datenum(1:di:end);
         rdat.time_cmp = rdat.datenum(1:di_cmp:end);
      else  % gusT
         di_cmp        = round(length(rdat.time)/length(rdat.compass)); 
         rdat.time_cmp = rdat.time(1:di_cmp:end);
      end

   % sample rate in Hz
   sr        = round(1/(median(diff(rdat.time)*3600*24)));
   sr_cmp    = round(1/(median(diff(rdat.time_cmp)*3600*24)));


   % for normal time vector
   Nf    = dt*sr;     % number of data points in 10 min window
   J{1}  = 1:length(rdat.time);
   I     = split_fragments(J, Nf, 0);  

   % for compass time vector
   Nfc   = dt*sr_cmp;
   J{1}  = 1:length(rdat.time_cmp);
   Ic    = split_fragments(J, Nfc, 0);  


%_____________________Pitot______________________
   % which variable carries the Pitot signal

      % find pitot data W or WP
       if isfield(rdat, 'W')
         dV1 = abs(nanmean(rdat.W)-2.02);
         dV2 = abs(nanmean(rdat.WP)-2.02);
         if dV1>dV2
            W  = rdat.W;
         else
            W  = rdat.WP;
         end
       else  
         dV1 = abs(nanmean(rdat.W2)-2.02);
         dV2 = abs(nanmean(rdat.W3)-2.02);
         if dV1>dV2
            W  = rdat.W2;
         else
            W  = rdat.W3;
         end
      end




%_____________________initialize quantities______________________
   Praw.time  = nan(1,length(I));
   Praw.P     = nan(1,length(I));
   Praw.W     = nan(1,length(I));
   Praw.cmp   = nan(1,length(I));
   if is_chipod
      Praw.T1    = nan(1,length(I));
      Praw.T2    = nan(1,length(I));
      Praw.vT1   = nan(1,length(I));
      Praw.vT2   = nan(1,length(I));
   else  % gusT
      Praw.T    = nan(1,length(I));
      Praw.vT   = nan(1,length(I));
   end

%_____________________loop through dt intervals______________________

for i=1:length(I)

   %-------------normal bulk stuff-----------------
      Praw.time(i)  = nanmean( rdat.time(I{i}) ); 
      Praw.P(i)     = nanmean( rdat.P(I{i}) ); 

   %----------------temperature--------

   if is_chipod
      Praw.T1(i)    = nanmean( rdat.T1(I{i}) ); 
      Praw.T2(i)    = nanmean( rdat.T2(I{i}) ); 
      Praw.vT1(i)   = nanvar( rdat.T1(I{i}) ); 
      Praw.vT2(i)   = nanvar( rdat.T2(I{i}) ); 
   else
      Praw.T(i)    = nanmean( rdat.T(I{i}) ); 
      Praw.vT(i)   = nanvar( rdat.T(I{i}) ); 
   end

   %---------------------Pitot voltage--------------
      wtmp        = W(I{i}) ;
      %  removing negative outlieres
      wtmp( wtmp<(nanmean(wtmp)-2*nanstd(wtmp)) ) = nan;
      Praw.W(i)   = nanmean( wtmp ); 


   %--------------------compass----------------------
      if is_chipod
         tmp = rdat.CMP(Ic{i})/10/180*pi; % transfor into rad
      else %gusT
         tmp = rdat.compass(Ic{i})/180*pi; % transfor into rad
      end
      tmp         = exp(1i * tmp);  % convert to complex plain
      Praw.cmp(i) = angle(nanmean(tmp))/pi*180;

end

%_____________________save data______________________
      sdir = [basedir 'proc' filesep 'Praw' filesep];
      [~, ~, ~] = mkdir(sdir);
      is1 = strfind(rfid,'_');  % find under score
      is2 = strfind(rfid,'.');  % find dot
      savestamp   = [rfid((is1):(is2)) 'mat'];
      
      % in cases where no 'raw_' preceeds the date for the raw filenames,
      % need different string for savestamp
      if isempty(is1)
          savestamp   = ['_' rfid(1:(is2)) 'mat'];
      end
      
      save([sdir filesep  'Praw' savestamp], 'Praw');
