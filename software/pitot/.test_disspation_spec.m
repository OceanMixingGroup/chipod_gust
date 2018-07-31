addpath(genpath('~/arbeit/matlab_tbx/plotting_matlab/'));

clear all;
close all;

N_spec_points  = 300;

alpha = 11/15;


eps_in  = 1e-6;
nu   = 2.3e-6;
eta_in  = (nu^3/eps_in)^.25;

dk     = 1e-2;
k      =  1e-2:dk:1e3;
ii_mes = find( k>=1e1 & k<=5e1 );
ii_mes = ii_mes(1:round(length(ii_mes)/N_spec_points):end); % downsample


D   = 2*nu*alpha*eps_in^(2/3)*k.^(1/3).*exp(-3/2*alpha*(eta_in*k).^(4/3));
noise =  10.^((rand(size(D)))-.6);  
noise =  1+(2*(rand(size(D)))-1);  
D_pure = D;
D = D.*noise;

   tic
      [eps, eta, epss]  = dissipation_spec(k(ii_mes), D(ii_mes));
      epss
   toc
   tic
      ii_mes2        =  1:length(ii_mes);
      [eps, eta, epss]  = dissipation_spec(k(ii_mes(ii_mes2)), D(ii_mes(ii_mes2)), eta);
   toc
epss./eps_in
  

%_____________________test noise______________________

%     test.time = 1:1000;
%     test.noise =  .3*rand(size(test.time));  
%     test.spd   = .1;

%     [E, f] = gappy_psd( test.spd+test.noise , 512, 100, 10);
%     k = 2*pi/test.spd*f;
%     D = 2*nu*k.^2.*E*test.spd/2/pi; 
%     [eps, eta, epss]  = dissipation_spec(k, D);
%     eps



figure
   loglog(k,D)
   hold all
   loglog(k(ii_mes),D(ii_mes), 'Linewidth',2)
   loglog(k,D_pure,'k', 'Linewidth',2)
