function [fig ] = plot_basic_temp(T, unit, dtind , tl, vis) 
%% [fig ] = plot_basic_temp(T, [unit], [dt] , [tl], [vis]) 
%     
%     Basic plotting routine for temp.mat data
%
%     INPUT 
%        T        :  temp.amt structure
%        unit     :  string containing unit
%        dt       :  step width for time array (default = 1)
%        tl       :  time limites
%        vis      :  shall the figure be visible (default = 'on')
%
%


%_____________________default values______________________
if nargin < 5
   vis = 'on';
end
if nargin < 4
   tl = T.time([1 end]);
end
if nargin < 4
   dtind = 1;
end
if nargin < 2
   unit = '';
end



   %tl, T, unit

   col = get(groot,'DefaultAxesColorOrder');
   

    fig = figure('Color',[1 1 1],'visible',vis,'Paperunits','centimeters',...
            'Papersize',[30 30],'PaperPosition',[0 -1 30 30]);


            [ax, ~] = create_axes(fig, 6, 1, 0);
            
            a = 1;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.depth(1:dtind:end), 'color', [0  0 0 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['P [dbar]'], 7);
               t = text_corner(ax(1), ['unit ' unit], -2);
               
            a = 2;
            if isfield(T, 'T1') 
               pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T1(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
               pj = 2; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T2(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                  t1 = text_corner(ax(a), ['T1 [^\circ C]'], 1);
                  t1.Color = [col(1,:)];
                  t2 = text_corner(ax(a), ['T2 [^\circ C]'], 3);
                  t2.Color = [col(2,:)];
            else
               pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T(1:dtind:end), 'color', [0  0 0 1], 'Linewidth', 1);
                  xlim(ax(a), tl);
                  t = text_corner(ax(a), ['T [^\circ C]'], 1);
            end
               xlim(ax(a), tl);
               
            
            a = 3;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.AX(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
            pj = 2; p(pj) = plot(ax(a), T.time(1:dtind:end), T.AY(1:dtind:end) + 4, ...
                                 'color', [col(pj,:) 1], 'Linewidth', 1);
            pj = 3; p(pj) = plot(ax(a), T.time(1:dtind:end), T.AZ(1:dtind:end) + 9.81 - 4, ...
                                 'color', [col(pj,:) 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t1 = text_corner(ax(a), ['AX [m s^{-2}]'], 1);
               t1.Color = [col(1,:)];
               t2 = text_corner(ax(a), ['AY (+4) [m s^{-2}]'], 2);
               t2.Color = [col(2,:)];
               t3 = text_corner(ax(a), ['AZ (+9.81 - 4) [m s^{-2}]'], 3);
               t3.Color = [col(3,:)];
               plot(ax(a), xlim(ax(a)), [1 1]*0, '--', 'color', col(1,:)*0.5);
               plot(ax(a), xlim(ax(a)), [1 1]*4, '--', 'color', col(2,:)*0.5);
               plot(ax(a), xlim(ax(a)), [1 1]*-4, '--', 'color', col(3,:)*0.5);
               ylim(ax(a), [-1, 1]*8);
            
            a = 4;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.cmp(1:dtind:end), '.', 'color', [0  0 0 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['compass [deg]'], 7);

            a = 5;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.W(1:dtind:end), 'color', [0  0 0 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['pitot [volt]'], 7);

            a = 6;
            if isfield(T, 'T1') 
              if isfield(T, 'varT1p')              
                pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.varT1p(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                pj = 2; p(pj) = plot(ax(a), T.time(1:dtind:end), T.varT2p(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                  xlim(ax(a), tl);
                  t1 = text_corner(ax(a), ['varT1P [K^2/s^2]'], 1);  
                  t1.Color = [col(1,:)];
                  t2 = text_corner(ax(a), ['varT2P [K^2/s^2]'], 3);  
                  t2.Color = [col(2,:)];
                  set(ax(a), 'yscale', 'log')
              else
               pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T1Pt(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
               pj = 2; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T2Pt(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                  xlim(ax(a), tl);
                  t1 = text_corner(ax(a), ['T1P [K/s]'], 1);  
                  t1.Color = [col(1,:)];
                  t2 = text_corner(ax(a), ['T2P [K/s]'], 3);  
                  t2.Color = [col(2,:)];
               Tpcombined = [T.T1Pt(1:dtind:end), T.T2Pt(1:dtind:end)];
               ylim(ax(a), 4*[-1, 1]*max(abs([prctile(Tpcombined, 2), prctile(Tpcombined, 98)])));
              end
            else  % gust
              if isfield(T, 'varTp')              
                  pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.varTp(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                     xlim(ax(a), tl);
                     t1 = text_corner(ax(a), ['varTP [K^2/s^2]'], 1);  
                     t1.Color = [col(1,:)];
                     set(ax(a), 'yscale', 'log')
              else
                  pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.TPt(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                     xlim(ax(a), tl);
                     t1 = text_corner(ax(a), ['TP [K/s]'], 1);  
                     t1.Color = [col(1,:)];
                   ylim(ax(a), 4*[-1, 1] * max(abs([prctile(T.TPt, 2), prctile(T.TPt, 98)])));
              end
            end
               
            linkaxes(ax, 'x');   

            abc='abcdefghijklmnopqrst';
            for a = 1:(size(ax,1)*size(ax,2))
               text_corner(ax(a), abc(a), 9);
               %set(ax(a), 'Xtick', ceil(tl(1)):round(diff(tl)/5):floor(tl(2)));
            end
            
            datetick(ax(a),  'keeplimits');
            if diff(tl)<2
               xlabel(ax(a), [datestr(floor(mean(tl)), 'dd mmm yyyy')]);
            elseif diff(tl)<35
               xlabel(ax(a), [datestr(floor(mean(tl)), 'mmm yyyy')]);
            else
               xlabel(ax(a), [datestr(floor(mean(tl)), 'yyyy')]);
            end
            
