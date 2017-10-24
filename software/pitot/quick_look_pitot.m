function [fig, P] = quick_look_pitot(data, W, vis)
%% [fig, P] = quick_look_pitot(data, W, [visible])
%
%  This function does a quick Pitot calibration based on the calibration coefs
%  in W including epsilon and spectrogram
%
% INTPUT
%  data      : standart output of quick_look, or quick_calibrate
%  W         : calibration structure out of header_p.mat
%  visible   : figure visible 'on' or ['off'];
%
% OUTPUT
%  P.time    : time-vector for speed
%  P.spd     : speed vector
%  P.U       : complex velocity vector
%  P.time_eps : time vector for epsilon
%  P.f        : frequency vector
%  P.eps_f    : normalized spectrogram(f,t) (unit epsilon)
%  P.eps      : epsilon(t)
%   
%   created by: 
%        Johannes Becherer
%        Tue Oct 24 09:45:26 PDT 2017

if nargin < 3
   vis = 'on';
end


 % calibrate Pitot tube
    P.time = data.time;

    [P.spd, ~, ~] = pitot_calibrate(data.W, data.T, data.P, W);

     spd_vel = P.spd; 
     spd_vel(P.spd<0) = 0; % remove negative speeds before add direction
     [P.U] = pitot_add_direction(P.time, spd_vel, data.time_cmp, data.cmp);

%_____________________calculate eps______________________

    [spec, f, T] = fast_spectrogram(P.time, P.spd, 1/(24*24), 1/(24*24)/2); 
    [ eps, var_eps, eps_f] = icscaling_velocity(T, f, spec, P.time, P.spd, [1/50 1/20]);

   P.time_eps = T;
   P.f        = f;
   P.eps_f    = eps_f;
   P.eps      = eps;  

%_____________________generate figure______________________

    fig = figure('Color',[1 1 1],'visible',vis,'Paperunits','centimeters',...
            'Papersize',[30 20],'PaperPosition',[0 0 30 20])
    

      [ax, ~] = create_axes(fig, 6, 1, 0);
         squeeze_axes(ax , .95, .88);
         shift_axes(ax, 0, .1);
         squeeze_axes(ax(6) , 1, 2);
         shift_axes(ax(6), 0, -.14);

      col = get(groot,'DefaultAxesColorOrder');
      xl = data.time([1 end]);
      
      a=1;
       if isfield(data, 'T')
          pj = 1;
         plot(ax(a), data.time, data.T, 'color', [col(pj,:) .5], 'Linewidth', 1);
         plot(ax(a), data.time, movmean(data.T, 60000), 'color', [col(pj,:) .5], 'Linewidth', 2);
      else
        pj = 1; 
         plot(ax(a), data.time, data.T1, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time, movmean(data.T1, 60000), 'color', [col(pj,:) .5], 'Linewidth', 2);
        pj = 2; 
         plot(ax(a), data.time, data.T2, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time, movmean(data.T2, 60000), 'color', [col(pj,:) .5], 'Linewidth', 2);
         legend(p, 'T1', 'T2');
      end
         xlim(ax(a), xl);
         t = text_corner(ax(a), ['temperature [deg C]'], 1);
         t = text_corner(ax(a), [datestr(median(data.time), 'dd mmm yyyy')], -2);
         
         
      
      a=2;
         plot(ax(a), P.time, P.spd, 'color', [col(pj,:) .5], 'Linewidth', 1);
         plot(ax(a), P.time, movmean(P.spd, 60000), 'color', [col(pj,:) .5], 'Linewidth', 2);
            xlim(ax(a), xl);
            t = text_corner(ax(a), ['speed [m/s]'], 1);

      a=3;
        pj = 1; 
         plot(ax(a), data.time_cmp, data.cmp, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time_cmp, movmean(data.cmp, round(60000/25)), 'color', [col(pj,:) .5], 'Linewidth', 2);
        pj = 2; 
         plot(ax(a), data.time_cmp, data.pitch, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time_cmp, movmean(data.pitch, round(60000/25)), 'color', [col(pj,:) .5], 'Linewidth', 2);
        pj = 3; 
         plot(ax(a), data.time_cmp, data.roll, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time_cmp, movmean(data.roll, round(60000/25)), 'color', [col(pj,:) .5], 'Linewidth', 2);
            xlim(ax(a), xl);
            t = text_corner(ax(a), ['orientation [deg]'], 1);
            legend(p, 'compass', 'pitch', 'roll');
         

      a=4;
        pj = 1; 
         plot(ax(a), data.time, data.a_vel_x, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time, movmean(data.a_vel_x, 60000), 'color', [col(pj,:) .5], 'Linewidth', 2);
        pj = 2; 
         plot(ax(a), data.time, data.a_vel_y, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time, movmean(data.a_vel_y, 60000), 'color', [col(pj,:) .5], 'Linewidth', 2);
        pj = 3; 
         plot(ax(a), data.time, data.a_vel_z, 'color', [col(pj,:) .5], 'Linewidth', 1);
         p(pj) = plot(ax(a), data.time, movmean(data.a_vel_z, 60000), 'color', [col(pj,:) .5], 'Linewidth', 2);
            xlim(ax(a), xl);
            ylim(ax(a), [-2 2]);
            t = text_corner(ax(a), ['vel_{acc} [m/s]'], 1);
            legend(p, 'vel_x', 'vel_y', 'vel_z');
         
         
      a=5;
      pj=1;
         plot(ax(a), P.time_eps, P.eps, 'color', [col(pj,:)], 'Linewidth', 2);
            xlim(ax(a), xl);
            t = text_corner(ax(a), ['\epsilon [m2/s3]'], 1);
            set(ax(a), 'Yscale', 'log');
            ylim(ax(a), [1e-10 1e-2]);

      a=6;
         pcolor(ax(a), P.time_eps, P.f, log10(P.eps_f))
            shading(ax(a), 'flat');
            caxis(ax(a), [-8 -2]);
            colormap(ax(a), 'jet')
            cb = colorbar('peer', ax(a));
					axpos = get(ax(a), 'Position');
					cb.Position =  [axpos(1)+axpos(3)+.01 axpos(2) .02 axpos(4) ];
				set(ax(a), 'Yscale', 'log');
               px = xl'; py = [1/50 1/20];
               patch([px px(2) px(1)], [py(1) py(1) py(2) py(2)], [.3 .3 .3], ...
                     'facealpha', .3, 'edgecolor', [.5 0 0], 'Linewidth', 1, 'parent', ax(a));
            t = text_corner(ax(a), ['\epsilon [m2/s3]'], 1);
            t.BackgroundColor = [1 1 1 .7];
            ylabel(ax(a), 'f [Hz]')
            xlim(ax(a), xl);
            
            
            abc='abcdefghijklmnopqrst';
            for a = 1:(size(ax,1)*size(ax,2))
               tabc = text_corner(ax(a), abc(a), 7);
               tabc.BackgroundColor = [1 1 1 .5];
               set(ax(a), 'box', 'on', 'TickDir', 'out');
            end
            

            datetick(ax(a), 'keeplimits');
            linkaxes(ax, 'x');
            

               
            
				
                  


         



         

