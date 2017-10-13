function [y1] =  clever_interp( x, y, x1 )
%% [y1] =  clever_interp( x, y, x1 ) 
%  
%  This function interpts the y(x) -> y1(x1) in a clever way:
%     1) taking care of nans
%     2) if  dt(x) < dt(x1) a corresponding low pass filter is applied
%
%  Note x and x1 are expected to be close to equidistant
%
%  
%   created by: 
%        Johannes Becherer
%        Tue Oct 10 13:33:31 PDT 2017


dt    = nanmean( diff(x) );
dt1   = nanmean( diff(x1));

if dt1 > 2*dt
   ww =  round(dt1/dt);
   y  =  movmean( y, ww);
end

if size(x,2) == size(y,1)
    y = y'; 
end

ii = ~(isnan(x) | isnan(y));

y1 = interp1( x(ii), y(ii), x1);

