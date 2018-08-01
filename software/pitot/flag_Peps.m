function [Peps]  = flag_Peps(Peps)
%% [Peps]  = flag_Peps(Peps)
%
%     This function contructs the flaggin masks for 
%     pitot epsilon structure
%     
%   created by: 
%        Johannes Becherer
%        Tue Jul 31 15:41:46 PDT 2018



   dt = median(diff(Peps.time));
   Dt = 20/3600/24;  %20sec
   Nlp = round(Dt/dt);

   Peps.spd_lp = movmean(Peps.spd, Nlp);
   Peps.spd_hp = Peps.spd - Peps.spd_lp;
   
   % mask for local maxima
   Peps.mask_max   = nan(size(Peps.time));
   Peps.mask_max(find(Peps.spd_hp>0)) = 1;
   Peps.mask_max   = movmean(Peps.mask_max, 3);

   % mask wave amplitude vs mean current
   Peps.mask_wave      = nan(size(Peps.time));
   Peps.wave_amplitude = sqrt(2)*movstd(Peps.spd, Nlp);
   Peps.mask_wave(find(log10(abs(Peps.wave_amplitude./Peps.spd_lp))<0)) = 1;
   Peps.mask_wave   = movmean(Peps.mask_wave, 2);


   % mask flow speed mim 10 cm/sec
   Peps.mask_flow      = nan(size(Peps.time));
   Peps.mask_flow(find(Peps.spd_lp>.1)) = 1;
   Peps.mask_flow   = movmean(Peps.mask_flow, 2);


   Peps.mask_tot = Peps.mask_max.*Peps.mask_wave.*Peps.mask_flow;

   Peps.eps_nomask = Peps.eps;
   Peps.eps    =  Peps.eps.*Peps.mask_tot;
end
