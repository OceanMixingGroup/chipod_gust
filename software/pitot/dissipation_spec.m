function [eps, eta, epss]  = dissipation_spec(k_in, D_in, nu, eta_in)
%% [eps, eta, epss]  = dissipation_spec(k, D, eta)
%
%     This function calculates epsilon based on the dissipation spectra 
%     (gradient spectra)  
%     D(k) = k^2*E(k)   
%     where E(k) is the Energy spectra of the velocity
%
%     INPUT:
%        f     : f - vector
%        D     : dissipation spectrum D(f)
%        nu    :  molecular disspation
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
if nargin < 3
   nu = 1.3e-6;
end

% iteration
Nmax_iteration = 10;
epss   = nan(1,Nmax_iteration);
eta =nan;
eps = nan;

D2eps = 15*nu;  % 15 because 1D longitudinal

Neta = 3; % how many kolmogorov length should be the kmax cut off

% assuming that k_in is equidistant
dk_in = median(diff(k_in));
kmin_in  =  k_in(1);
kmax_in  =  k_in(end);
ii_in    =  1:length(k_in);

% integral over input k_range
%eps_in  = sum(D_in)*dk_in; 
%eps_in  = sum(D_in)*dk_in; 
eps_in  = D2eps*int_eps( k_in, D_in);
if nargin<4
   eta_in = (nu^3/(eps_in*2))^.25;
end

   % ignor high waves numbers since they might be contaminated
   if kmax_in > 1/(Neta*eta_in)
      kmax_in = 1/(Neta*eta_in);
     % eps_in  = sum(D_in(k_in<=kmax_in))*dk_in; 
      ii_in = find(k_in<=kmax_in);

      if length(ii_in) > 2
         eps_in  = D2eps*int_eps( k_in(ii_in), D_in(ii_in));
      else  % no enough spectral points to work with
         return;
      end
   end





epss(1) = eps_in;

[G1_na,kks]=nasmyth_G1(1000,20);

   for i = 1:Nmax_iteration;
      
      [D_na, k_na, eta] = nasmyth2Dk( kks, G1_na, epss(i), nu );

         if k_in(ii_in(end)) > 1/(Neta*eta)
            kmax_in = 1/(Neta*eta);
            ii_in = find(k_in<=kmax_in);
            if length(ii_in) > 2
               eps_in  = D2eps*int_eps( k_in(ii_in), D_in(ii_in));
            else  % no enough spectral points to work with
               return
            end
         end

      ii_na_part  =  find( k_na>= kmin_in & k_na<= kmax_in);
      eps_na_part = D2eps*int_eps( k_na(ii_na_part), D_na(ii_na_part)); 
      eps_na      = D2eps*int_eps( k_na, D_na); 



       epss(i+1)  = eps_na*(eps_in/eps_na_part);


      % break for early convergence
      if abs(epss(i+1)/epss(i)-1) < .01 % 99% convergence 
         break;
      end
   end

%_____________________Output______________________
ii_nan = find(~isnan(epss));
eps = epss(ii_nan(end));
eta = (nu^3/eps)^.25;

   % only for testing purpose
   if 0%  rand(1)>.99
        loglog(k_in, D_in);
        hold all;
        loglog( k_in(k_in<=kmax_in), D_in(k_in<=kmax_in), 'g'); 
        loglog(k_na, D_na);
        loglog(k_na*2*pi, D_na/2/pi, '--');
        loglog([1 1]*1/eta,[1e-5 1e-3], 'r--' )
        loglog([1 1]*1/(5*eta),[1e-5 1e-3], 'r--' )
        t = text_corner(gca, [num2str(epss([1 ii_nan(end)]))], 1);
        keyboard;
        hold off;
   end


end

function [Ix] =  int_eps(k,D)
   [~, whichdim] = max(size(k));
   Ix = trapz( k, D, whichdim);
end
