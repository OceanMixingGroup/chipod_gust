%     this is meant to study individual raw files 
%
%   created by: 
%        Johannes Becherer
%        Tue Nov 22 10:21:01 PST 2016

clear all;
close all;

%_____________________flags______________________
do_adcp           =  0;  % do you have ADCP data to compare with
header_offset     =  0;  % use offset data from pitot header
get_v0_hist       =  0;  % detremine v0 from histogram
get_v0_fit        =  0;  % detremine v0 from adcp fit
do_vel_plot       =  0;  % make a big velocity comparision figure
do_eps            =  0;  % calculate disspation rates from Pitot

ifid               = 1;  % which raw files you would like to look at


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines

%_____________________load ADCP (or other velocity data)______________________
if do_adcp
   load ../input/vel_m.mat;
   A.U   = vel_m.u + 1i*vel_m.v;
   A.spd = abs(A.U); 
   A.time = vel_m.time;
end

%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   
   rfid = fids{ifid};


%_____________________calibrate data______________________
   [data, head] = quick_look([basedir '/raw/'], rfid);


   %_____________________teperature______________________
   if isfield(data, 'T1')
      % if chipod pick temperature sensor closed to pitot
      data.T = data.T1; 
   else
      disp(['is gusT']);
   end


   %--------------------calibrate pitot----------------------
   %load([ddir '/calib/head_p.mat']);
   load([basedir '/calib/header_p.mat']);

   if header_offset
      % off-sets from header
        data.cal.V0 = W.V0;
        data.cal.T0 = W.T0
        data.cal.P0 = W.P0
   else
      % maual pick off sets
        data.cal.V0 = 0;
        data.cal.T0 = 27;
        %data.cal.T0 = nanmedian(data.T);
        data.cal.P0 = 0;
        %data.cal.P0 = nanmedian(data.P);
   end

    % slopes from header
        data.cal.Vs = 1/W.Pd(2);
        data.cal.Ts = W.T(2);
        data.cal.Ps = W.Ps(2);

   % do standart calibration
   if header_offset
        [data.spd, data.Pdym, data.V_cal] = pitot_calibrate(data.W, data.T, data.P,...
                 data.cal.V0, data.cal.T0, data.cal.P0, data.cal.Vs, data.cal.Ts, data.cal.Ps);
   end


%_____________________determine V0______________________

   %---------------------via histogram----------------------

   if get_v0_hist & ~header_offset  
       % do pre calibration
         [~, ~, pre.V_cal] = pitot_calibrate(data.W, data.T, data.P,...
                    data.cal.V0, data.cal.T0, data.cal.P0, data.cal.Vs, data.cal.Ts, data.cal.Ps);

       % determine V0
         V0_hist = pitot_v0_hist(pre.V_cal, data.cal, 1)
             data.cal.V0 = V0_hist;

       % repeat calibration with new V0
        [data.spd, data.Pdym, data.V_cal] = pitot_calibrate(data.W, data.T, data.P,...
                 data.cal.V0, data.cal.T0, data.cal.P0, data.cal.Vs, data.cal.Ts, data.cal.Ps);

       % add direction ot Pitot
       [data.U] = pitot_add_direction(data.time, data.spd, data.time_cmp, data.cmp);
   end


    %---------------------via ADCP fit----------------------

    if get_v0_fit & do_adcp  & ~header_offset  
       % do pre calibration
         [~, ~, pre.V_cal] = pitot_calibrate(data.W, data.T, data.P,...
                    data.cal.V0, data.cal.T0, data.cal.P0, data.cal.Vs, data.cal.Ts, data.cal.Ps);

       % determine V0
          [V0_adcp] = fit_pitot_v0( A.time, A.spd, data.time, pre.V_cal, data.cal.Vs, 1)
          data.cal.V0 = V0_adcp;

       % repeat calibration with new V0
        [data.spd, data.Pdym, data.V_cal] = pitot_calibrate(data.W, data.T, data.P,...
                 data.cal.V0, data.cal.T0, data.cal.P0, data.cal.Vs, data.cal.Ts, data.cal.Ps);

       % add direction ot Pitot
       [data.U] = pitot_add_direction(data.time, data.spd, data.time_cmp, data.cmp);

    end



%_____________________cal epsilon______________________
   if do_eps
      tic
      [eps, data1, M] = chi_cal_pitot_eps(data, W);
      toc


      ii = (abs(eps.slope+5/3)<.5);

      figure
         plot( eps.time, eps.eps, 'Linewidth', 1);
         hold all;
         ww=10;
         plot( moving_median(eps.time, ww, ww/5), moving_average(eps.eps, ww, ww/5), 'Linewidth', 2);
          plot( moving_median(eps.time(ii), ww, ww/5), moving_average(eps.eps(ii), ww, ww/5), 'k','Linewidth', 2);

   end


%_____________________velocity comparision______________________

 if do_vel_plot & do_adcp
      a_time = A.time;
      a_filt = diff(a_time(1:2))*3600*24/600 % 10 min low pass
      a_U    = qbutter( A.U, a_filt);
      a_L    = 'ADCP';
      p_time = data.time;
      p_filt = diff(p_time(1:2))*3600*24/600 % 10 min low pass
      p_U    = qbutter(data.U, p_filt);
      p_L    = 'Pitot';
      [fig] =  compare_velocity_timeseries(a_time, a_U, a_L, p_time, p_U, p_L);
      %print(gcf,'./chi_vs_ADCP.png','-dpng','-r200','-painters')
  end
   



