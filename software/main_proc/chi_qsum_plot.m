function []  = chi_qsum_plot(basedir, vis)
%%  chi_qsum_plot(basedir, vis)
%     
%     This function generates a summary figure


   load([basedir 'proc' filesep 'qsum.mat']);

   tl = S.time([1 end]);


   col = get(groot,'DefaultAxesColorOrder');
   

    fig = figure('Color',[1 1 1],'visible',vis,'Paperunits','centimeters',...
            'Papersize',[30 30],'PaperPosition',[0 -1 30 30])


            [ax, ~] = create_axes(fig, 5, 1, 0);
            
            a = 1;
            pj = 1; p(pj) = plot(ax(a), S.time, S.depth, 'color', [0  0 0 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['P [dbar]'], 7);
               
            a = 2;
            if isfield(S, 'T1') 
               pj = 1; p(pj) = plot(ax(a), S.time, S.T1, 'color', [col(pj,:) 1], 'Linewidth', 2);
               pj = 2; p(pj) = plot(ax(a), S.time, S.T2, 'color', [col(pj,:) 1], 'Linewidth', 2);
                  t1 = text_corner(ax(a), ['T1 [^\circ C]'], 1);
                  t1.Color = [col(1,:)];
                  t2 = text_corner(ax(a), ['T2 [^\circ C]'], 3);
                  t2.Color = [col(2,:)];
            else
               pj = 1; p(pj) = plot(ax(a), S.time, S.T, 'color', [0  0 0 1], 'Linewidth', 2);
                  xlim(ax(a), tl);
                  t = text_corner(ax(a), ['T [^\circ C]'], 1);
            end
               xlim(ax(a), tl);
               
            
            a = 3;
            pj = 1; p(pj) = plot(ax(a), S.time, S.AX, 'color', [col(pj,:) 1], 'Linewidth', 2);
            pj = 2; p(pj) = plot(ax(a), S.time, S.AY, 'color', [col(pj,:) 1], 'Linewidth', 2);
            pj = 3; p(pj) = plot(ax(a), S.time, S.AZ, 'color', [col(pj,:) 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t1 = text_corner(ax(a), ['AX [m s^{-2}]'], 1);
               t1.Color = [col(1,:)];
               t2 = text_corner(ax(a), ['AY [m s^{-2}]'], 2);
               t2.Color = [col(2,:)];
               t3 = text_corner(ax(a), ['AZ [m s^{-2}]'], 3);
               t3.Color = [col(3,:)];
            
            a = 4;
            pj = 1; p(pj) = plot(ax(a), S.time, S.cmp, 'color', [0  0 0 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['compass [deg]'], 7);

            a = 5;
            pj = 1; p(pj) = plot(ax(a), S.time, S.W, 'color', [0  0 0 1], 'Linewidth', 2);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['pitot [volt]'], 7);

            linkaxes(ax, 'x');   

            abc='abcdefghijklmnopqrst';
            for a = 1:(size(ax,1)*size(ax,2))
               text_corner(ax(a), abc(a), 9);
               set(ax(a), 'Xtick', ceil(tl(1)):round(diff(tl)/5):floor(tl(2)));
            end
            


            datetick(ax(a), 'keepticks',  'keeplimits');
            
            unit    = chi_get_unit_name(basedir); % get unit name
            t = text_corner(ax(1), ['unit ' unit], -2);
            
            
      

%_____________________save pic______________________

print(gcf,[basedir 'pics' filesep 'qsum.png' ],'-dpng','-r200','-painters')
print(gcf,[basedir 'pics' filesep 'qsum.pdf' ],'-dpdf','-painters')



