function [Tz_eff, Gamma, Kt]  =  Tz_bbl_eff( eps, chi, mab, T, S, Depth)


   kappa =  .4;
   g     = 9.81;
   alpha = sw_alpha(S,T,Depth, 'temp');

   Kt       =  eps.^(1/3) *(kappa*mab)^(4/3); 

   Tz_eff   =  2^-.5 * chi.^.5 .*Kt.^-.5; 

   Gamma    = g *(.5)^.5 * alpha.*chi.^.5 .*eps.^(-5/6) *(kappa*mab)^(2/3);


end
