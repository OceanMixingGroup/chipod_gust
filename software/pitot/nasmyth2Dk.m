function [D, k, eta] = nasmyth2Dk( kks, G1, eps, nu )
%%   [D, k, eta] = nasmyth2Dk( kks, G1, eps, nu )
%     
%        this function dimensionalizes the spectrum from
%        [G1, kks] = nasmyth_G1(1000,20)
%
%        D     :   gradient spectrum k^2Phi_11  ->  eps = 15*nu*int D dk
%                    [s^-2 /(rad/m)]
%
%        k     :  radian wave number [rad/m]
%        eta   :  Kolmogorov length  
%
%   created by: 
%        Johannes Becherer
%        Tue Jan 29 10:24:48 PST 2019
   
if nargin < 4
   nu = 1.3e-6
end

      eta   = (nu^3/eps)^.25;
      k  =  kks/eta*2*pi;  
      D  =  G1*eta^(-2)*(eps*nu^5)^.25/(2*pi);
end 
