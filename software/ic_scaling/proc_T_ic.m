function [Tic] =  proc_T_ic( Ttime, T, Tztime, Tz, N2, Utime, U, specLength, DeltaTime, f_range, saveSpec )
%% [Tic] =  proc_T_ic( Ttime, T, Tztime, Tz, N2, Utime, U, specLength, DeltaTime, f_range, saveSpec )
%
%  This function make an IC batchlor scaling for a given temperature time series
%
%  INPUT
%     Ttime       :  time vector of temperature
%     T           :  Temperature vector
%     Tztime      :  time vector of stratification
%     Tz          :  vertical temrepature gradient
%     N2          :  N2 (same time step as Tz)
%     Utime       :  time vector for velocity
%     U           :  velocity vector (can be complex)
%     specLength  :  length of the spectrum (default 1/(24*6) = 10 min )
%     DeltaTime   :  time step of scaling (default = specLength)
%     f_range     :  frequncy range for the scaling (default [.02 .05])
%     saveSpec    :  boolean (if 1 save the spectrogram data, 0 don't [default])
%
%  OUTPUT
%     Tic.time    :  time vector for fit (dt = DeltaTime)
%     Tic.T       :  mean temerpature 
%     Tic.Tz      :  mean vertical tempreature gradient
%     Tic.N2      :  mean vertical stratification
%     Tic.U       :  mean velocity vector
%     Tic.chi     :  tempreature variance decay 
%     Tic.eps     :  dissipatoin of turbulent kinetic energy
%     
%     if savespec
%     Tic.f       :  frequncy vector of spectrogram
%     Tic.spec    :  spectrogram of temperature 
%     Tic.eps_f   :  normalized spectrogram [m^2s^-3]
%     Tic.chi_f   :  normalized spectrogram [K^2s^-3]
%
%
%   created by: 
%        Johannes Becherer
%        Fri Oct 13 10:49:12 PDT 2017

%---------------------set default values----------------------
if nargin < 8
   specLength  =  1/(24*6);
end
if nargin < 9
   DeltaTime   =  specLength;
end
if nargin < 10
   f_range   =  [.02 .05];
end
if nargin < 11
   saveSpec  =  0;
end
gamma = .2;
c_tau = .4; 

%_take care of dimensions_
if size(Ttime,1) > size(Ttime,2)
    Ttime = Ttime';
end
if size(T,1) > size(T,2)
    T = T';
end

%_____________________calculate spectrogram of temperature______________________
   [spec, f, Tic.time] = fast_spectrogram( Ttime, T, specLength, DeltaTime);

%_____________________interpolate on spectrum time step______________________
   Tic.T   = clever_interp( Ttime, T, Tic.time );
   Tic.U   = clever_interp( Utime, U, Tic.time );
   Tic.Tz  = clever_interp( Tztime, Tz, Tic.time );
   Tic.N2  = clever_interp( Tztime, N2, Tic.time );

%_____________________normalize spectrogram______________________

   %eps_f    =  (Tic.N2.*(Tic.Tz.^-2)/2/gamma/c_tau).^1.5 * (2*pi)^-2.*abs(Tic.U).^-1 ...
   %                  .* spec.^1.5 .* f'.^2.5;
   %chi_f    =  (Tic.N2.*(Tic.Tz.^-2)/2/gamma).^.5*c_tau^-1.5 * (2*pi)^-2.*abs(Tic.U).^-1 ...
   %                  .* spec.^1.5 .* f'.^2.5;
       eps_f    =  (Tic.N2.*(Tic.Tz.^-2)/2/gamma/c_tau).^1.5* (2*pi).* abs(Tic.U).^-1 ...
                     .* spec.^1.5 .* f'.^2.5;
   chi_f    =  (Tic.N2.*(Tic.Tz.^-2)/2/gamma).^.5*c_tau^-1.5 * (2*pi).*abs(Tic.U).^-1 ...
                     .* spec.^1.5 .* f'.^2.5;

%_____________________average eps and chi in f_range______________________

   if f_range(1)< f(1)
      f_range(1) = f(1);
   end
   iif      = find( f>=f_range(1) & f<=f_range(2) );
   Tic.eps  = mean( eps_f( iif,: ),1);
   Tic.chi  = mean( chi_f( iif,: ),1);


%_____________________extra output______________________
Tic.f_range =  f_range;

if saveSpec
   Tic.f       = f;
   Tic.spec    = spec;
   Tic.eps_f   = eps_f;
   Tic.chi_f   = chi_f;
end




