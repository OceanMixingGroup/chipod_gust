function [chi] = chi_chi_proc_ic(T, S, Tz, Tw, f_range)
%% [chi] = chi_chi_proc_ic(T, S, Tz, [Tw], [f_range])
%
%   This is the main chi processing routine
%
%   INPUT
%     T.time   :  time vector of T
%     T.T      :  T
%     [T.S]    :  S (default set to 35 psu) 
%     T.depth  :  depth
%     S.time   :  time vector of spd 
%     S.spd    :  speed
%     Tz.time  :  time vector of Tz and N2   
%     Tz.Tz    :  Tz
%     Tz.N2    :  N2
%     Tw       :  length of spectral window (default 600sec)
%     f_range  :  fitting range (default 1/50...1/20)
%
%   OUTPUT
%     chi.chi  :  chi data
%     chi.eps  :  dissipation rates
%     chi.time :  time vector
%     chi.spd  :  average speed used for the calculation
%     chi.T    :  temperature used to calculate viscosity
%     chi.S    :  salinity used to calculate viscosity
%     chi.depth:  depth used to calculate viscosity
%     chi.Tz   :  temperature gradietn used to do the calculations 
%     chi.N2   :  N2 used to do the calculations
%     chi.mask :  mask used  
%
%   created by: 
%        Johannes Becherer
%        Mon Jan 30 17:14:53 PST 2017
  

%_____________________default values______________________
if nargin < 4
   Tw    = 600;            % length of the window in sec
end
if ~isfield(T, 'S') 
   T.S = ones(size([T.T]))*35;
end

% fitting range 
if nargin < 5
   f_range  =  [1/50 1/20];
end
if diff(f_range) < 0
   f_range = f_range([2 1]);
end
fstart = f_range(1); fstop = f_range(2);  % integration range

    

%_____________________bring spd and T to same time grid______________________
   if length(S.time) ~= length(T.time)
      S.spd  = interp1( S.time, S.spd, T.time, 'nearest');
      S.time = T.time;
   end


%_____________________cut time series into pieces for spectra______________________
   dt    = nanmedian(diff(T.time))*3600*24;
   Nf    = round(Tw/(dt));   % Nf is the length of the Tw fragment
   Nf_2  = round(Nf/2);      % half a segment

   ii_s  =  [(1+Nf_2):Nf:(length(T.time)-Nf_2)]; % index array for chi-subsample


   J{1}  = 1:length(T.time); 
   I     = split_fragments(J, Nf, 0); % split in Tw intrevals with no overlap  
   Ni    = length(I);   % total number of Tw fragments  


%_____________________construct final time vector______________________
   chi.time    = T.time(ii_s);
   chi.spd     = clever_interp(S.time, S.spd, chi.time);
   chi.dTdz    = clever_interp( Tz.time, Tz.Tz, chi.time);
   chi.N2      = clever_interp( Tz.time, Tz.N2, chi.time);
   chi.T       = clever_interp( T.time, T.T, chi.time);
   chi.S       = clever_interp( T.time, T.S, chi.time);
   chi.depth   = clever_interp( T.time, T.depth, chi.time);


   %----------calculate viscosity and diffusivity-------------
      nu   = nan(Ni, 1);
      tdif = nan(Ni, 1);
      for i = 1:Ni
         nu(i)    =  sw_visc( chi.S(i), chi.T(i), chi.depth(i));
         tdif(i)  =  sw_tdif( chi.S(i), chi.T(i), chi.depth(i));
      end
   %---------------------spectral constants----------------------
         samplerate = Nf/Tw;
         nfft       = floor(Nf/4)*2; % use two windows on entire time length making sure nfft is even

%_____________________masking data______________________
   % avg.fspd(ik) >= 0.04  && dTdz(i)>min_dTdz 
   chi.mask = ones(size(chi.time));

   chi.mask(isnan(chi.dTdz))        = 0;
   chi.mask(isnan(chi.N2))          = 0;
   chi.mask(chi.N2<0)               = 0;



%_____________________major processing ______________________
   %---------------------initialize----------------------
   chi.chi     = nan(1,Ni); 
   chi.eps     = nan(1,Ni);

   %-------loop through all 1 sec seqments---------------
   for i = 1:Ni

      if chi.mask(i)  % check if calculation should be executed 

          [t_power, freq] = fast_psd( T.T(I{i}) ,nfft, samplerate);

          % calculate psd of dT/dt 
          tp_power   = t_power.*(2*pi*freq).^2;
                               

          % fit the spectrum
         [chi_tmp, eps_tmp , k,spec, k_kraich, spec_kraich, stats] =...
                    get_chipod_chi_ic( freq,fstart,fstop, tp_power, chi.spd(i), nu(i), tdif(i), chi.dTdz(i), chi.N2(i));


          chi.chi(i) = chi_tmp(1);
          chi.eps(i) = eps_tmp(1);
                     

      else % masked out 
          chi.chi(i)  =  nan;
          chi.eps(i)  =  nan;
      end
            

   end

