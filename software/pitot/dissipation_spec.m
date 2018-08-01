function [eps, eta, epss]  = dissipation_spec(k_in, D_in, eta_in)
%% [eps, eta, epss]  = dissipation_spec(k, D, [kmax])
%
%     This function calculates epsilon based on the dissipation spectra 
%     (gradient spectra)  
%     D(k) = 2*nu*k^2*E(k)   
%     where E(k) is the Energy spectra of the velocity
%
%     INPUT:
%        k     : k - vector
%        D     : dissipation spectrum D(k)
%        eta_in : maximum wave number for the fit
%                (default is 1e4 which is very slow T ~ kmax
%                but provides sufficient estimates 
%                for eps<1e-2 )
%                
%                for quick stratigy start with kmax = 1e-2
%                 and if eta < 5/kmax 
%                        kmax = 5/eta
%                        redo
%                 
%     OUTPUT:
%        eps   :  epsilon disspation rates
%        eta   :  koresponding kolmogorov length
%        epss  :  all iteration steps for diagnose purpose
%
%
%   created by: 
%        Johannes Becherer
%        Tue Jul 10 16:56:14 PDT 2018




% parameters
nu = 2.3e-6;
alpha = 11/15; % kolmogorov constant

% assuming that k_in is equidistant
dk_in = median(diff(k_in));
kmin_in  =  k_in(1);
kmax_in  =  k_in(end);

% integral over input k_range
eps_in  = sum(D_in)*dk_in; 
if nargin<3
   eta_in = (nu^3/(eps_in*2))^.25;
end

% ignor high waves numbers since they might be contaminated
if kmax_in > 5/eta_in
   kmax_in = 5/eta_in;
   eps_in  = sum(D_in(k_in<=kmax_in))*dk_in; 
end


      

   


% for model based on Pao 

   % should be at least 10x smaller than kmin_input;
   % 3x smaller than Delta k input (to get good resolution)
   % 100x smaller than 1/eta
   kmin  = min( [kmin_in*1e-1, dk_in/3, 1/eta_in*1e-2] ); 

   % should be at least 10x bigger than kmax_in;
   % and 5 bigger than 1/eta
   kmax  = max( [kmax_in*10, 5/eta_in]); 
   dk    = kmin;
   k = kmin:dk:kmax;
   ii_in =  find( k>=(kmin_in -dk_in*.5) & k<=(kmax_in+dk_in*.5));




% iteration
Nmax_iteration = 10;
epss   = nan(1,Nmax_iteration);
epss(1) = eps_in;
   for i = 1:Nmax_iteration;
      [D_tmp, eta]   =  pao_spectrum(k, epss(i));
      intDkmin = 3/2*nu*alpha*epss(i)^(2/3)*kmin.^(4/3);
      epss(i+1)  = ( sum(D_tmp)*dk + intDkmin )/( sum(D_tmp(ii_in))*dk )*epss(1);

      % break for early convergence
      if abs(epss(i+1)/epss(i)-1) < .01 % 99% convergence 
         break;
      end
   end



%_____________________Output______________________
ii_nan = find(~isnan(epss));
eps = epss(ii_nan(end));
eta = (nu^3/eps)^.25;

