function [D_off] = find_cmp_offset(basedir)
%%      [D_off] = find_cmp_offset(basedir)
%     This functnio calculates the compass offset based on
%     temp.mat and vel_m.mat 
%     does not need Pitot velocity
%
%   created by: 
%        Johannes Becherer
%        Thu Mar 21 17:23:04 PDT 2019

D_off = 0;

Tfid    =[basedir '/proc/temp.mat'];
vel_fid = [basedir '/input/vel_m.mat'];
if ~exist(Tfid)
   disp('I could not calculate compass offset');
   disp('you have to generate temp.mat first!');
   return ;
end
if ~exist(vel_fid)
   disp('I could not calculate compass offset');
   disp('you have to generate vel_m.mat first!');
   return ;
end

load( Tfid);
load( vel_fid);

[TL] =   whoAmI_timeLimits(basedir);

T.Uadcp   =  clever_interp(vel_m.time, vel_m.U, T.time);
[T.Um]   =  pitot_add_direction(T.time, abs(T.Uadcp), T.time_cmp, T.cmp);
tl(1) = max( [T.time(1) vel_m.time(1) TL.master(1)]);
tl(2) = min( [T.time(end) vel_m.time(end) TL.master(2)]);

   
   ww = 60*60*20*median(diff(T.time))*3600*24;
   D_off = nanmedian(angle(movmean(T.Um, ww)) - angle( movmean(T.Uadcp, ww)))*180/pi;

   disp('The calculated direction off set between ADCP and Chipod/GusT is');
   disp([num2str( D_off ) ' deg']);

