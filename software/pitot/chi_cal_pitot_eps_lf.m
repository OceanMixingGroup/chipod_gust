function [eps, varargout] = chi_cal_pitot_eps_lf(data, W)
%% [eps, [data, M]] = chi_cal_pitot_eps_lf(data, W) 
%     
%     This function calculates dissipation rates based on Pitot speeds for low frequencies < wave band
%
%     INPUT
%        data     :  calibrated chipod/ gusT data
%        data.W   :  must be contained Raw-voltage of Pitot
%        W        :  pitot header should contain (V0, T0, P0, Pd, T, Ps) [like W] 
%
%     OUTPUT
%        eps      :  product structure 
%                       contains : time, eps, slope, l_k, m_spd, M_spd, M_turn, M_tot
%        data     :  like input structure with additional velocity information 
%                       spd, spd_e, U, ...
%
%
%   created by: 
%        Johannes Becherer
%        Tue Jan 31 17:04:23 PST 2017




%_____________________for chipods______________________
if isfield(data, 'T1')
   data.T = data.T1;
end

%_____________________calibrate speed______________________
   [data.spd, data.Pdym, data.V_cal] = pitot_calibrate(data.W, data.T, data.P, W);

   % get directional information for U
   data.U  = pitot_add_direction( data.time, data.spd, data.time_cmp, data.cmp);

%_____________________filter speed______________________

   %---------------------define masks----------------------------------
      
         % no masks for now

   %---------------------apply filter----------------------
   %data.spd_e           = data.spd;
   %data.spd_e(M.spd)    = nan;
   %data.spd_e(M.turn)   = nan;
   
   

%_____________________cal epsilon______________________

   % frequency range
   f_range = [1/200 1/20]; % lower than wave band

%_____________________cut time series into 10 min pieces for spectra______________________
   Tw    = 600;            % length of the window in sec
   dt    = diff(data.time(1:2))*3600*24;
   Nf    = round(Tw/(dt));   % Nf is the length of the Tw fragment
   J{1}  = 1:length(data.time); 
   I     = split_fragments(J, Nf, 0); % split in Tw intrevals with no overlap  
   Ni    = length(I);   % total number of Tw fragments  
%_____________________construct final time vector______________________
   eps.time = nan(1,Ni);
   eps.spd  = nan(size(eps.time));
   for i = 1:Ni
      eps.time(i) = nanmean( data.time( I{i}) );
       % average speed
       eps.spd(i) = nanmean(data.spd(I{i}));
   end

   %---------------------initialize----------------------
   eps.eps      = nan(1,length(I));
   eps.slope    = nan(1,length(I));
   eps.l_k      = nan(1,length(I));

   %---------------------fit data----------------------
   for i = 1:Ni
      warning('off')
      [eps.eps(i), eps.slope(i), ~, eps.l_k(i)] = fit_longitudinal_kolmogorov(data.time(I{i}),  data.spd(I{i}), f_range, 0);
   end


   
   
   %_____________________additional output______________________
   varargout{1} = data;
