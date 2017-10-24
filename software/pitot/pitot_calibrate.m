function [spd, Pdym, Vcal] = pitot_calibrate(Vraw, T, P, W, T0, P0, Vs, Ts, Ps)
%% [spd, Pdym, Vcal] = pitot_calibrate(Vraw, T, P, W)
%   Converts raw pitot voltages into usful data
%
%   OUTPUT
%      spd     : Speed [m/s]
%      Pdym    : dynamic pressure [Pa]
%      Vcal    : calibrated pitot Voltage [V]
%
%   INPUT
%      Vraw    : raw pitot Voltage [V]
%      T       : temperature time sieries corresponding to Vraw (should have same length as Vraw)
%      P       : pressure time sieries corresponding to Vraw (can have same length as Vraw or be a scalar)
%      W.T0      : Temperature at corresponding V0 [C]
%      W.P0      : Pressure at corresponding V0 [psi] 
%      W.V0      : minimum voltage [V]
%      W.Ts      : temperature calibration 
%      W.Ps      : static Pressure calibration 
%      W.Pd      : dynamic pressure calibration 
%
%   created by: 
%        Johannes Becherer
%        Wed Nov 16 11:02:39 PST 2016
%%

if nargin < 5
    T0 = W.T0;
    P0 = W.P0;
    Vs = 1/W.Pd(2);
    Ts = W.T(2);
    Ps = W.Ps(2);
    V0 = W.V0;
else
    V0 = W;
end
    


% If P is a scalar (for mooring deployments)
if(length(P)<2) 
   P = ones(size(Vraw))*P;
end
if(length(V0)<2) 
   V0 = ones(size(Vraw))*V0;
end
if(length(T)<2) 
   T = ones(size(Vraw))*T;
end

Sp = Vs; % slope for the dynamic pressure
    if Sp <0
        Sp = 1/Sp;
    end
ST = Ts; % slope for the temperature
SP = Ps; % slope for the static pressure


Vcal = Vraw - (T-T0)*ST - (P-P0)*SP - V0;
Pdym = Vcal*Sp;
spd    = sign(Pdym).*sqrt(2/1025*abs(Pdym)); 
