clear all;
close all;

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines

%____________________ load data ______________________

load ../input/vel_m.mat;
   A.time = vel_m.time;
   A.U    = vel_m.u + 1i*vel_m.v;

load ../proc/pitot.mat;
   P.time = pitot.time;
   P.U    = pitot.U;

%____________________ bring on the same time step ______________
adt = diff(A.time(1:2));
pdt = diff(P.time(1:2));

a_filt = adt/pdt;
p_filt = pdt/adt;

if a_filt > 1
   a_filt = 1;
end

if p_filt > 1
   p_filt = 1;
end


% integer step

astep = round(1/a_filt);
pstep = round(1/p_filt);

% butter filter
      a_U    = qbutter( A.U, a_filt);
      p_U    = qbutter( P.U, p_filt);


%_____________________velocity comparision______________________

      a_L    = 'ADCP';
      p_L    = 'Pitot';
      [fig] =  compare_velocity_timeseries(A.time(1:astep:end), A.U(1:astep:end), a_L, P.time(1:pstep:end), P.U(1:pstep:end), p_L);
      print(gcf,'../pics/Pitot_vs_ADCP.png','-dpng','-r200','-painters')
