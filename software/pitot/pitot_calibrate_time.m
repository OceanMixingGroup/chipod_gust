function [spd, Pdym, Vcal] = pitot_calibrate_time(time, Vraw, T, P,  W)
%% [spd, Pdym, Vcal] = pitot_calibrate_time(time, Vraw, T, P, W)
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
%        Wed Aug 15 20:53:54 PDT 2018

    T0 = W.T0;
    P0 = W.P0;
    Vs = 1/W.Pd(2);
    Ts = W.T(2);
    Ps = W.Ps(2);
    V0 = W.V0;

if(length(V0)<2) 
   V0 = ones(size(Vraw))*V0;
else
   V0 = interp1( W.time, W.V0, time, 'nearest', 'extrap');
end

[spd, Pdym, Vcal] = pitot_calibrate(Vraw, T, P, V0, T0, P0, Vs, Ts, Ps);

