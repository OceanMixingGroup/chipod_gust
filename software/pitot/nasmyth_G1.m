function[sp,k1]= nasmyth_G1(n,shrink);
% [sp,k1]= nasmyth_G1(n,shrink);
%
%  Adapted from org Jim's fitting routine but using Lueck adapted 
%  (by me,.. not sure if there is anything official) fitting
%   for G1 instead of G2 
%
%   this is meant for longitudinal spectra from pitot instead of
%   transfers spectra from chameleon.
%
%        Johannes Becherer
%        Mon Jan 28 13:08:37 PST 2019
% 
%
% % old description from nasmyth.m 
% program to evaluate the spectral coefficients according to
% Oakey's spline fit to the G2 form of the empirical Nasmyth spectrum 
% sp gives the spectral values at frequencies k1.  n equally spaced spectral
% values are  given up to a maximum of  (k/ks)=2pi/shrink

% where:
%       n is the number of spectral estimates you want back (it is optional,
%       th)
%       k1: k/ks, where k is the cyclic waveno., ks is Kolm. waveno [rad/m]
%       [k1]=rad^-1
%       sp: non-dimensional spectrum form G2(k/ks)
%
%               these 2 variables can now be scaled to dimensional
%               quantities ug
%              eps, nu
%       ref: Moum, May 1996

if nargin<1
  n=1000;
  shrink=1;
end

ii=(1/n/shrink):(1/n/shrink):1/shrink;
k1=2*pi*ii'; %[rad^-1]

sp  = 6.02*k1.^(1/3)./(1+(27.5*k1).^3.7);

return

	
