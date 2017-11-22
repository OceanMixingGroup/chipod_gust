%% this script does all the basic calibration excluding chi-processing and save 1s infomation
%     in ./proc/temp.mat
%  
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:01:26 PDT 2017


do_parallel = 1;     % use paralelle computing 
do_raw_proc = 1;     % do the raw data processing 
do_plot     = 0;     % generate a over view plot 

time_range = [datenum(2000, 1, 1, 0, 0, 0) ...
              datenum(2060, 1, 1, 0, 0, 0)];

dtind = 600; % every 10 minutes, assuming 1 second estimates

% parameters for spectra of T1, T2
nfft = 30*600; % length of segment (in seconds)
nbandsmooth = 1; % number of frequency bands to average over

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location


if do_raw_proc
   
   % generate temp.mat
   generate_temp( basedir, do_parallel, time_range)

   %____________________create motion.mat file______________________

   load([basedir 'proc' filesep 'temp.mat']);

   motion.chiangle = atan2(T.a_vel_y, T.a_vel_x) * 180/pi;
   motion.time = T.time;
   motion.cmp = T.cmp;
   motion.a_vel_x = T.a_vel_x;
   motion.a_vel_y = T.a_vel_y;
   motion.a_vel_z = T.a_vel_z;

   save([basedir 'proc' filesep 'motion.mat'], 'motion', '-v7.3');


end

%_____________________make multiple summary plots______________________
if do_plot

  if ~exist('T', 'var')
      load([basedir 'proc' filesep 'temp.mat']);
  end

   tind(1) = find_approx(T.time, time_range(1), 1);
   tind(2) = find_approx(T.time, time_range(2), 1);
   tl = T.time(tind);
   unit    = chi_get_unit_name(basedir); % get unit name

   fig1 = plot_basic_temp(T, unit, dtind, tl);


	%_____________________save pic______________________
	print(fig1,[basedir 'pics' filesep 'temp.png' ],'-dpng','-r200','-painters')
   %	savefig(gcf,[basedir 'pics' filesep 'temp.fig' ])

        %% displacement histogram

        CreateFigure; hold on;
        histogram(sqrt(2)*movstd(T.a_dis_x(tind(1):tind(2)), 60), ...
                  'displaystyle', 'stairs', 'linewidth', 1.5)
        histogram(sqrt(2)*movstd(T.a_dis_y(tind(1):tind(2)), 60), ...
                  'displaystyle', 'stairs', 'linewidth', 1.5)
        histogram(sqrt(2)*movstd(T.a_dis_z(tind(1):tind(2)), 60), ...
                  'displaystyle', 'stairs', 'linewidth', 1.5)
        legend('x', 'y', 'z');
        xlabel('Displacement = sqrt(2) x std(integrated accel_{x,y,z}) over 1 minute')
        print(gcf,[basedir 'pics' filesep 'disp.png' ],'-dpng','-r200','-painters')

        %% spectra of T1 and T2
        fs = 1/(diff(T.time(1:2))*86400);

        trange = [find_approx(T.time, time_range(1)):find_approx(T.time, time_range(2))];

        tic; disp('Calculating PSD')
        if isfield(T, 'T1')
           [p1,f1] = fast_psd(T.T1(trange), nfft, fs);
           [p2,f2] = fast_psd(T.T2(trange), nfft, fs);

           % frequency band smoothing
           f = moving_average(f1, nbandsmooth, nbandsmooth);
           p1a = moving_average(p1, nbandsmooth, nbandsmooth);
           p2a = moving_average(p2, nbandsmooth, nbandsmooth);
        else
           [p1,f1] = fast_psd(T.T(trange), nfft, fs);
           f = moving_average(f1, nbandsmooth, nbandsmooth);
           p1a = moving_average(p1, nbandsmooth, nbandsmooth);
        end
        toc;

        tscale = 1; %in seconds
        CreateFigure;
        loglog(1./f/tscale, p1a); hold on;
        if isfield(T, 'T1')
           loglog(1./f/tscale, p2a);
           legend('T1', 'T2');
        else
           legend('T');
        end
        set(gca, 'XDir', 'reverse')
        xlabel('Period (s)')
        ylabel('PSD');
        title([unit ' | nfft = ' num2str(nfft/60) 's | nbandsmooth = ' num2str(nbandsmooth)])

        print(gcf,[basedir 'pics' filesep 'temp-spectra.png' ],'-dpng','-r200','-painters')

        %% histograms comparing chipod body motion to background flow
        CreateFigure;

        dt = diff(T.time(1:2))*86400;
        trange = [find_approx(T.time, time_range(1)):find_approx(T.time, time_range(2))];

        histogram(abs(T.a_vel_x(trange)), 'normalization' ,'pdf', ...
                  'displaystyle', 'stairs', 'linewidth', 1.5, ...
                  'DisplayName', '|v_x|');
        hold on;
        histogram(hypot(T.a_vel_x(trange), T.a_vel_y(trange)), ...
                  'normalization' ,'pdf', 'displaystyle', 'stairs', ...
                  'linewidth', 1.5, 'displayname', '|v_x + v_y|');
        histogram(abs(T.a_vel_z(trange)), 'normalization' ,'pdf', ...
                  'displaystyle', 'stairs', 'linewidth', 1.5, ...
                  'displayname', '|v_z|');

        titlestr = [num2str(round(dt)) 's avg for v_*'];

        try
            load ../input/vel_m.mat

            trange2 = [find_approx(vel_m.time, time_range(1)):find_approx(vel_m.time, time_range(2))];
            histogram(vel_m.spd(trange2), 'normalization' ,'pdf', ...
                      'displaystyle', 'stairs', 'linewidth', ...
                      1.5, 'displayname', 'adcp');
            dt2 = diff(vel_m.time(1:2))*86400/60;
            titlestr = [titlestr ' | ' ...
                        num2str(round(dt2)) 'minutes \Deltat for background'];
        catch ME
            disp('could not load vel_m.mat');
            disp('skipping');
        end

        try
            load ../input/vel_p.mat

            trange2 = [find_approx(vel_p.time, time_range(1)):find_approx(vel_p.time, time_range(2))];
            histogram(vel_p.spd(trange2), 'normalization' ,'pdf', ...
                      'displaystyle', 'stairs', 'linewidth', ...
                      1.5, 'displayname', 'pitot');
        catch ME
            disp('could not load vel_p.mat');
            disp('skipping');
        end

        ylim([0 10]);
        xlim([0 1.5]);
        xlabel('m/s')
        ylabel('PDF')
        title(titlestr)
        legend('-dynamiclegend');

        print(gcf,['../pics/velocity-histogram.png'],'-dpng','-r200','-painters')
end
