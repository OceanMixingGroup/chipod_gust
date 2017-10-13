function [spec, f, T] = fast_spectrogram(t, x, SWidth, TSteps)
%% [spec, f, T] = fast_spectrogram(t, x, SWidth, TSteps)
%        This function generates a spectrogram PS of the vector x
%         x       : 1 dimensional data vector
%         t       : coresponding time (equidistant and in days)
%         SWidth  : spectrum width in days
%         TSteps  : Time step for the spectrogram in days.

%   %test
%   SWidth = 1000;
%   TSteps = 499;
%   t = 1:.5:3000;
%   x = sin(t/200*2*pi);

   dt  = nanmean(diff(t)) ; % time step in days
   Nww = floor(SWidth/dt);
     
   % % set spectrum parameters
   %NFFTmax = Nww;       % set largest FFT-segment 
   NFFTmax = 2^floor(log2(Nww));       % set largest FFT-segment 
   fs       =  1/(dt*3600*24);


   IT = 1:Nww;
   cnt = 1;
   DT = round(TSteps/dt);
   ii = IT;
   while ii(end)<length(t)
       T(cnt) = mean(t(ii));
       %[f , psd(:,cnt), ~, ~]  = iow_fancypsd(x(ii), NFFTmax, Ffac, Fmin, fs, p, pPlot);
       [spec(:, cnt), f] = fast_psd( x(ii) , NFFTmax, fs);
       ii = [ii] + DT;
       cnt = cnt+1;
   end
