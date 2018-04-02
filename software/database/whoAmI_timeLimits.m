function [TL] =   whoAmI_timeLimits(basedir)
%%    [TL] =   whoAmI_timeLimits(basedir)
%
%     This function extracts all usefull time limits from whoAmI.m
%
%     INPUT
%        basedir     :  unit directory
%
%     OUTPUT   
%        TL          :  structure with all time limits
%
%
%   created by: 
%        Johannes Becherer
%        Tue Mar 27 12:40:09 PDT 2018
%

% check fi whoAMI exist?
fname =  [basedir '/mfiles/whoAmI.m'];

tl_default   = [datenum(2000, 1,1) datenum(2050, 1,1)];  % default tl

%_____________________initialize______________________

TL.master   =  tl_default;
TL.T1       =  tl_default;
TL.Tp1      =  tl_default;
TL.T2       =  tl_default;
TL.Tp2      =  tl_default;
TL.P        =  tl_default;
TL.pitot    =  tl_default;
TL.acc      =  tl_default;
TL.cmp      =  tl_default;


%_____________________get from whoAmI______________________
if exist(fname)

  addpath(fname);
  DB = whoAmI;

  TL.master(1)       =  convert_timeStr2datenum( DB.instruments.start, tl_default(1));
  TL.master(2)       =  convert_timeStr2datenum( DB.instruments.stop, tl_default(2));

  tl_default = TL.master;


  %_____________________set all sub flags to master______________________
   TL.T1       =  TL.master;
   TL.Tp1      =  TL.master;
   TL.T2       =  TL.master;
   TL.Tp2      =  TL.master;
   TL.P        =  TL.master;
   TL.pitot    =  TL.master;
   TL.acc      =  TL.master;
   TL.cmp      =  TL.master;




  % chipod
   if isfield(DB, 'chipods')

      TL.T1(2)    =  convert_timeStr2datenum( DB.chipods.T1_stop, tl_default(2));
      TL.T2(2)    =  convert_timeStr2datenum( DB.chipods.T2_stop, tl_default(2));
      TL.Tp1(2)   =  convert_timeStr2datenum( DB.chipods.Tp1_stop, tl_default(2));
      TL.Tp2(2)   =  convert_timeStr2datenum( DB.chipods.Tp2_stop, tl_default(2));

      TL.P(2)     =  convert_timeStr2datenum( DB.chipods.P_stop, tl_default(2));
      TL.pitot(2) =  convert_timeStr2datenum( DB.chipods.pitot_stop, tl_default(2));
      TL.acc(2)   =  convert_timeStr2datenum( DB.chipods.acc_stop, tl_default(2));
      TL.cmp(2)   =  convert_timeStr2datenum( DB.chipods.cmp_stop, tl_default(2));

   else  % gust

      TL.T1(2) =  convert_timeStr2datenum( DB.gust.T_stop, tl_default(2));
      TL.Tp1(2) =  convert_timeStr2datenum( DB.gust.Tp_stop, tl_default(2));

      TL.P(2)     =  convert_timeStr2datenum( DB.gust.P_stop, tl_default(2));
      TL.pitot(2) =  convert_timeStr2datenum( DB.gust.pitot_stop, tl_default(2));
      TL.acc(2)   =  convert_timeStr2datenum( DB.gust.acc_stop, tl_default(2));
      TL.cmp(2)   =  convert_timeStr2datenum( DB.gust.cmp_stop, tl_default(2));

   end


end

end


function [t] = convert_timeStr2datenum( timeStr, defaultTime )

   if ~(isempty(timeStr))
         try 
            t = datenum( timeStr);
         catch
            t = defaultTime;
         end
         return;
   end

   t = defaultTime;
   

end

