function [data] = chi_calibrate_chipod(rfid, head)
%%
%     
%        This function reads CHIPOD raw data and calibrates them acording
%        to the coefficients in head.
%
%        Input
%           rfid   :  path the the specific raw-file
%           head   :  corresponding header
%
%        Output
%           data   : data structure containing calibrated data
%
%   created by: 
%        Johannes Becherer
%        Fri Sep  2 15:53:26 PDT 2016


%_____________________read in raw data______________________

[rdat, headtest]  = raw_load_chipod(rfid);

%____________ check if raw_file head corresponds to given head___________

if ~isequal(head.coef, headtest.coef)
   disp(' ');
   disp('!!! WARNING !!!!!');
   disp('the raw file header is different from the given header');
   disp('the given header is used and the raw-file header ignored');
   disp(' ');
end

%_____________________calibrate data______________________


   % tempertuare
      chi.T1=calibrate_polynomial(rdat.T1,head.coef.T1);
      chi.T2=calibrate_polynomial(rdat.T2,head.coef.T2);
   % pressure
      chi.P=calibrate_polynomial(rdat.P,head.coef.P);
      chi.depth = (chi.P-14.7)/1.47;

   % time vector
      chi.time    = rdat.datenum(1:2:end);
      chi.time_tp  = rdat.datenum;
      %chi.time_cmp = chi.time(1:10:end);
      chi.time_cmp = chi.time(1:10:(round(length(chi.time)*.1)*10));
      
   % accelerometer
         g=9.81;
         chi.AX=g.*calibrate_polynomial(rdat.AX,head.coef.AX);
         chi.AY=g.*calibrate_polynomial(rdat.AY,head.coef.AY);
         chi.AZ=g.*calibrate_polynomial(rdat.AZ,head.coef.AZ);
          chi.AX=fillgap(chi.AX);
          chi.AY=fillgap(chi.AY);
          chi.AZ=fillgap(chi.AZ);
          [dis,vel]=integrate_acc(chi,head);

          chi.a_dis_x = dis.x;
          chi.a_dis_y = dis.y;
          chi.a_dis_z = dis.z;
          chi.a_vel_x = vel.x;
          chi.a_vel_y = vel.y;
          chi.a_vel_z = vel.z;
           
         chi.AXtilt=calibrate_tilt(rdat.AX,head.coef.AX);
         chi.AYtilt=calibrate_tilt(rdat.AY,head.coef.AY);
         chi.AZtilt=calibrate_tilt(rdat.AZ,head.coef.AZ);

   % compass
         chi.cmp = rdat.CMP/10+head.coef.CMP(1);
         % In case time_cmp is to long (sometimes its one index to long)
         chi.time_cmp = chi.time_cmp(1:length(chi.cmp));
      

   % DTdt
      rdat.T1P = rdat.T1P - nanmean(rdat.T1P);
      rdat.T2P = rdat.T2P - nanmean(rdat.T2P);
      chi.T1Pt = calibrate_tp( rdat.T1P, head.coef.T1P, rdat.T1, head.coef.T1, 100*ones(size(rdat.T1)) );
      chi.T2Pt = calibrate_tp( rdat.T2P, head.coef.T2P, rdat.T2, head.coef.T2, 100*ones(size(rdat.T1)) );


   %---------------------Pitot stuff----------------------

   % pitot_voltage
      % find pitot data W or WP
      chi.W = pitot_choose_W(rdat);

      [chi.spd, chi.Pdym, chi.V_cal] = pitot_calibrate(chi.W, chi.T1, chi.P, head.W);

      chi.U  = pitot_add_direction( chi.time, chi.spd, chi.time_cmp, chi.cmp);
      chi.u  = real(chi.U);
      chi.v  = imag(chi.U);
         



%---------------------return data----------------------
data = chi;


