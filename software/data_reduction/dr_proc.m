 function [chi] =  dr_proc(DR, head, vel)
  % clear all;
  % close all;

%  addpath(genpath('../../software/'));% include  path to preocessing routines
%  load ../../../../proc/DR/dr__14060100.mat;
%  load ../../../../proc/chi/chi_mi11_ic/chi_14060100.mat;
%     chi_r = chi;
%     clear chi;
%  load ../../../../input/vel_m.mat;
%  load ../../../../calib/header.mat;
%     vel.time = vel_m.time;
%     vel.spd  = vel_m.spd;

   


g = 9.81;

%_____________________calibrate basics______________________

   % time
   chi.time = DR.time;

   % temperature
   chi.T1   =  (DR.T1.^2+DR.vT1)*head.coef.T1(3)+ DR.T1*head.coef.T1(2) + head.coef.T1(1);
   chi.T2   =  (DR.T2.^2+DR.vT2)*head.coef.T2(3)+ DR.T2*head.coef.T2(2) + head.coef.T2(1);
 
   % pressure
   chi.P    =  DR.P*head.coef.P(2) + head.coef.P(1);
   chi.depth=  chi.P/1.47;

   % accelleration
   chi.AX    =  g*(DR.AX*head.coef.AX(2) + head.coef.AX(1));
   chi.AY    =  g*(DR.AY*head.coef.AY(2) + head.coef.AY(1));
   chi.AZ    =  g*(DR.AZ*head.coef.AZ(2) + head.coef.AZ(1));

   chi.vAX   =  DR.vAX*head.coef.AX(2)^2*g^2;
   chi.vAY   =  DR.vAY*head.coef.AY(2)^2*g^2;
   chi.vAZ   =  DR.vAZ*head.coef.AZ(2)^2*g^2;

   % chipod motion variance
   chi.vUax   =  DR.vUax*head.coef.AX(2)^2*g^2;
   chi.vUay   =  DR.vUay*head.coef.AY(2)^2*g^2;
   chi.vUaz   =  DR.vUaz*head.coef.AZ(2)^2*g^2;

   
   % calibration coeff for Tp
   Ctp1  = 2*DR.T1*head.coef.T1(3)+ head.coef.T1(2); 
   Ctp2  = 2*DR.T2*head.coef.T2(3)+ head.coef.T2(2); 

   % dTdz
   chi.Tz1  =  -DR.Tz1.*Ctp1/head.coef.P(2)*1.47;
   chi.Tz2  =  -DR.Tz2.*Ctp2/head.coef.P(2)*1.47;

   % N2
   alpha    = sw_alpha(35, mean(chi.T1), 30, 'temp');
   chi.N2_1 = g*chi.Tz1*alpha;
   chi.N2_2 = g*chi.Tz2*alpha;

%_____________________speed passed the sensor____________________

   % interpolate speed on DR time step
   chi.spd_b = interp1(vel.time, vel.spd, DR.time);
   % add motion variance
   chi.spd   =  sqrt(chi.spd_b.^2 + chi.vUax + chi.vUay + chi.vUaz);


   
%_____________________Turbulence calculations______________________


   % N2
         nu    =  sw_visc( 35, mean(chi.T1), 30);
         tdif  =  sw_tdif( 35, mean(chi.T1), 30);
   % Coeff spec
        Cpsd_1 =  (Ctp1/head.coef.T1P(2)).^2;
        Cpsd_2 =  (Ctp2/head.coef.T1P(2)).^2;


   % loop through 10min bits
   for i = 1:length(DR.f2)
         fstart = DR.f2{i}(1);
         fstop  = DR.f2{i}(end);

      % full spectrum
         % T1
         [chi_tmp, eps_tmp , k,spec, k_kraich, spec_kraich, stats] =...
                    get_chipod_chi_ic( DR.f1{i},fstart,fstop, DR.Pt1_1{i}*Cpsd_1(i) , chi.spd(i), nu, tdif, chi.Tz1(i), chi.N2_1(i));
            chi.chi1_1(i) = chi_tmp(1);
            chi.eps1_1(i) = eps_tmp(1);
         % T2
         [chi_tmp, eps_tmp , k,spec, k_kraich, spec_kraich, stats] =...
                    get_chipod_chi_ic( DR.f1{i},fstart,fstop, DR.Pt2_1{i}*Cpsd_2(i)  , chi.spd(i), nu, tdif, chi.Tz2(i), chi.N2_2(i));
            chi.chi2_1(i) = chi_tmp(1);
            chi.eps2_1(i) = eps_tmp(1);

      % cut spectrum
         % T1
         [chi_tmp, eps_tmp , k,spec, k_kraich, spec_kraich, stats] =...
                    get_chipod_chi_ic( DR.f2{i},fstart,fstop, DR.Pt1_2{i}*Cpsd_1(i) , chi.spd(i), nu, tdif, chi.Tz1(i), chi.N2_1(i));
            chi.chi1_2(i) = chi_tmp(1);
            chi.eps1_2(i) = eps_tmp(1);


      % cut spectrum
         % T1
         [chi_tmp, eps_tmp , k,spec, k_kraich, spec_kraich, stats] =...
                    get_chipod_chi_ic( DR.f3{i},fstart,fstop, DR.Pt1_3{i}*Cpsd_1(i) , chi.spd(i), nu, tdif, chi.Tz1(i), chi.N2_1(i));
            chi.chi1_3(i) = chi_tmp(1);
            chi.eps1_3(i) = eps_tmp(1);
         % T2
         [chi_tmp, eps_tmp , k,spec, k_kraich, spec_kraich, stats] =...
                    get_chipod_chi_ic( DR.f3{i},fstart,fstop, DR.Pt2_3{i}*Cpsd_2(i)  , chi.spd(i), nu, tdif, chi.Tz2(i), chi.N2_2(i));
            chi.chi2_3(i) = chi_tmp(1);
            chi.eps2_3(i) = eps_tmp(1);

      % slope 
         % T1
         chi.chi1_4(i) = ((DR.fit_Tp1(i)*Cpsd_1(i))/.4)^(3/2) * sqrt( chi.N2_1(i)/chi.Tz1(i)^2/2/.2 ) * (2*pi)^(-2)/chi.spd(i) ;
         chi.eps1_4(i) = ((DR.fit_Tp1(i)*Cpsd_1(i))/.4)^(3/2) * ( chi.N2_1(i)/chi.Tz1(i)^2/2/.2 )^(3/2)* (2*pi)^(-2)/chi.spd(i) ;

         chi.chi2_4(i) = ((DR.fit_Tp2(i)*Cpsd_2(i))/.4)^(3/2) * sqrt( chi.N2_2(i)/chi.Tz2(i)^2/2/.2 ) * (2*pi)^(-2)/chi.spd(i) ;
         chi.eps2_4(i) = ((DR.fit_Tp2(i)*Cpsd_2(i))/.4)^(3/2) * ( chi.N2_2(i)/chi.Tz2(i)^2/2/.2 )^(3/2)* (2*pi)^(-2)/chi.spd(i) ;
   end

