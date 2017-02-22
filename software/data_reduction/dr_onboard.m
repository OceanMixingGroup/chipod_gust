function  [DR] = dr_onboard(basedir, rfid, varargin)
%%    [DR] = dr_onboard(basedir, rfid, [dt])
%
%        This function emulates the processing on board of a 
%        chipod for data reduction. It returns a structure that
%        contains all the data that need to be send back in order
%        to garantee turbulence processing
%
%        INPUT
%           basedir  :  processing base directory
%           rfid     :  name of the raw file to process
%           dt       :  what average interval should be chosen in sec (default 600s)
%
%        OUTPUT
%           DR.time  :  time vector
%           DR.T1    :  mean temperture (sensor 1)
%           DR.T2    :  mean temperture (sensor 2)
%           DR.vT1   :  variance of temperture (sensor 1)
%           DR.vT2   :  variance of temperture (sensor 2)
%           DR.P     :  pressure
%           DR.Tz1   :  vertical temperature gradient
%           DR.Tz2   :  vertical temperature gradient
%           DR.AX    :  average acceleration in x direction
%           DR.AY    :  average acceleration in y direction
%           DR.AZ    :  average acceleration in z direction
%           DR.vAX   :  variance of acceleration in x direction
%           DR.vAY   :  variance of acceleration in y direction
%           DR.vAZ   :  variance of acceleration in z direction
%
%           DR.f1    :  frequency vector for the full spectrum
%           DR.Pt1_1 :  power spectral density for the derivative spectrum (sensor 1) 
%           DR.Pt2_1 :  power spectral density for the derivative spectrum (sensor 2) 
%
%           DR.f2    :  frequency vector for the cut spectrum
%           DR.Pt1_2 :  power spectral density for the derivative spectrum (sensor 1) 
%           DR.Pt2_2 :  power spectral density for the derivative spectrum (sensor 2) 
%
%           DR.f3    :  frequency vector for the bin averaged spectrum
%           DR.Pt1_3 :  power spectral density for the derivative spectrum (sensor 1) 
%           DR.Pt2_3 :  power spectral density for the derivative spectrum (sensor 2) 
%
%           DR.fit_Tp1 :  linear fit for psd (sensor 1) 
%           DR.fit_Tp1 :  linear fit for psd (sensor 2) 
%
%   created by: 
%        Johannes Becherer
%        Wed Feb  8 16:49:46 PST 2017

%  set averaging interval
if nargin < 3
   dt = 600;
else
   dt = varargin{1};
end


%_____________________read raw_data______________________
[rdat, ~] = raw_load_chipod([basedir '/raw/', rfid]);


%_____________________split data vector______________________
   % time vector
      rdat.time    = rdat.datenum(1:2:end);
      rdat.time_tp  = rdat.datenum;
      %rdat.time_cmp = rdat.datenum(1:10:end);

   % for normal time vector
   %Nf    = dt/median(diff(rdat.time)*3600*24);
   Nf    = dt/.02;
   J{1}  = 1:length(rdat.time);
   I     = split_fragments(J, Nf, 0);  

   % for Tp time vector
   %Nfp   = dt/median(diff(rdat.time_tp)*3600*24);
   Nfp   = dt/.01;
   J{1}  = 1:length(rdat.time_tp);
   Ip    = split_fragments(J, Nfp, 0);  

%_____________________spectral parametres______________________
   nfft       = Nfp/2; % use two windows on entire time length
   samplerate = 1/.01;
   f_min      = 1/200;
   f_max      = 1/20;

%_____________________Pitot______________________
   % which variable carries the Pitot signal
   W = rdat.W2; % only true for Pirata should be generalized

%_____________________initialize quantities______________________
   DR.time  = nan(1,length(I));
   DR.T1    = nan(1,length(I));
   DR.T2    = nan(1,length(I));
   DR.vT1   = nan(1,length(I));
   DR.vT2   = nan(1,length(I));
   DR.P     = nan(1,length(I));
   DR.Tz1   = nan(1,length(I));
   DR.Tz2   = nan(1,length(I));
   DR.AX    = nan(1,length(I));
   DR.AY    = nan(1,length(I));
   DR.AZ    = nan(1,length(I));
   DR.vAX   = nan(1,length(I));
   DR.vAY   = nan(1,length(I));
   DR.vAZ   = nan(1,length(I));
   DR.vUax   = nan(1,length(I));
   DR.vUay   = nan(1,length(I));
   DR.vUaz   = nan(1,length(I));
   DR.Wa     = nan(1,length(I));
   DR.Wa2    = nan(1,length(I));
   DR.Wm     = nan(1,length(I));

   DR.fit_Tp1 = nan(1,length(I));
   DR.fit_Tp2 = nan(1,length(I));

%_____________________integrate acceleration______________________
    % dummy head for integrate routine (This part should be integrated from the routine)
    dhead.sensor_index.AX = 1;
    dhead.samplerate = 50;
    [dis,vel]=integrate_acc(rdat,dhead);

%_____________________loop through dt intervals______________________

for i=1:length(I)

   %-------------normal bulk stuff-----------------
   DR.time(i)  = nanmean( rdat.time(I{i}) ); 
   DR.T1(i)    = nanmean( rdat.T1(I{i}) ); 
   DR.T2(i)    = nanmean( rdat.T2(I{i}) ); 
   DR.P(i)     = nanmean( rdat.P(I{i}) );
   DR.AX(i)    = nanmean( rdat.AX(I{i}) ); 
   DR.AY(i)    = nanmean( rdat.AY(I{i}) ); 
   DR.AZ(i)    = nanmean( rdat.AZ(I{i}) ); 

%_____________________Pitot______________________
   % which variable carries the Pitot signal

      % find pitot data W or WP
       if isfield(rdat, 'W')
         dV1 = abs(nanmean(rdat.W)-2.02);
         dV2 = abs(nanmean(rdat.WP)-2.02);
         if dV1>dV2
            W  = rdat.W;
         else
            W  = rdat.WP;
         end
       else  
         dV1 = abs(nanmean(rdat.W2)-2.02);
         dV2 = abs(nanmean(rdat.W3)-2.02);
         if dV1>dV2
            W  = rdat.W2;
         else
            W  = rdat.W3;
         end
      end
   wtmp        = W(I{i}) ;
   %  removing negative outlieres
   wtmp( wtmp<(nanmean(wtmp)-2*nanstd(wtmp)) ) = nan;

   DR.W(i)   = nanmean( wtmp ); 



   DR.vT1(i)   = nanvar( rdat.T1(I{i}) ); 
   DR.vT2(i)   = nanvar( rdat.T2(I{i}) ); 
   DR.vAX(i)   = nanvar( rdat.AX(I{i}) ); 
   DR.vAY(i)   = nanvar( rdat.AY(I{i}) ); 
   DR.vAZ(i)   = nanvar( rdat.AZ(I{i}) ); 

   %_____________________chipod velocity______________________
   DR.vUax(i)  = var(vel.x(I{i}));
   DR.vUay(i)  = var(vel.y(I{i}));
   DR.vUaz(i)  = var(vel.z(I{i}));

   %-----------temperature gradient---------------------
   %p = polyfit( rdat.P(I{i}), rdat.T1(I{i}), 1); % pressure has a drift problem
   p = polyfit( dis.z(I{i}), rdat.T1(I{i}), 1);
   DR.Tz1(i)  = p(1); 
   %p = polyfit( rdat.P(I{i}), rdat.T2(I{i}), 1);
   p = polyfit( dis.z(I{i}), rdat.T1(I{i}), 1);
   DR.Tz2(i)  = p(1); 

   %---------------------spectra----------------------

   % full (also integrate fast_psd in routine)   
   [DR.Pt1_1{i}, DR.f1{i}] = fast_psd( rdat.T1P(Ip{i}), nfft, samplerate);
   [DR.Pt2_1{i}, ~]        = fast_psd( rdat.T2P(Ip{i}), nfft, samplerate);
        
   % cut 
   iif = find( DR.f1{1}>=f_min & DR.f1{i}<= f_max);
   DR.f2{i}    = DR.f1{i}(iif);
   DR.Pt1_2{i} = DR.Pt1_1{i}(iif);
   DR.Pt2_2{i} = DR.Pt2_1{i}(iif);


   % bin average 
   Nf2 =  floor(length(DR.f2{i})/3);
   for j = 1:Nf2+1
      if j<=Nf2
         DR.f3{i}(j)    = mean( DR.f2{i}([1:3]+3*(j-1)) );
         DR.Pt1_3{i}(j) = mean( DR.Pt1_2{i}([1:3]+3*(j-1)) );
         DR.Pt2_3{i}(j) = mean( DR.Pt2_2{i}([1:3]+3*(j-1)) );
      else %combine rest of data points
         if length(DR.f2{i})/3 > Nf2 % only if there are data points left
            DR.f3{i}(j)    = mean( DR.f2{i}( (1+3*(j-1)):end) );
            DR.Pt1_3{i}(j) = mean( DR.Pt1_2{i}( (1+3*(j-1)):end) );
            DR.Pt2_3{i}(j) = mean( DR.Pt2_2{i}( (1+3*(j-1)):end) );
         end
      end
   end
    %  DR.f3    
    %  DR.Pt1_3 
    %  DR.Pt2_3 i
     %   iif = round(10.^((0:.2:log10(length(P100Hz.f)))));
     %   iif = unique(iif);
     %   [tmpf , ~, ~, ~]=databin( P100Hz.f(iif) , P100Hz.f , log10(P100Hz.f));



    % fit as f^(1/3) slope
    DR.fit_Tp1(i) = mean(DR.Pt1_2{i}./DR.f2{i}.^(1/3));
    DR.fit_Tp2(i) = mean(DR.Pt2_2{i}./DR.f2{i}.^(1/3));

end

