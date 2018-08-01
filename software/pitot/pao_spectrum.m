function [D, eta]   =  pao_spectrum(k, eps)
% [D, eta]   =  pao_spectrum(k, eps)
%   This function return the disspation spectrum
%   for a given wave number vector k and and
%   disspation rate of eps. Based on the Theoretical 
%   spectrum of Pao 65
%
%   created by: 
%        Johannes Becherer
%        Thu Jul 12 15:25:59 PDT 2018

% parameters
nu = 2.3e-6;
alpha = 11/15; % kolmogorov constant


eta  = (nu^3/eps)^.25;
D    = 2*nu*alpha*eps^(2/3)*k.^(1/3).*exp(-3/2*alpha*(eta*k).^(4/3));
