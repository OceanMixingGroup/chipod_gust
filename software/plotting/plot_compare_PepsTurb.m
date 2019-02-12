function [fig, fig1] = plot_compare_PepsTurb( basedir, whichTurb)
%%   [fig, fig1] = plot_compare_PepsTurb( basedir)
%
%  This function generate a basic diagnostic plot for the epsilon from Pitot tube
%
%  INPUT
%     basedir   :  instrument directory
%     whichTurb : default (mm1)
%
%   created by: 
%        Johannes Becherer
%        Thu Nov 30 14:54:08 PST 2017
if nargin < 2
   whichTurb = 'mm1';
end

   unit    = chi_get_unit_name(basedir); % get unit name

   load([basedir filesep 'proc' filesep 'pitot_eps.mat'])
   load([basedir filesep 'proc' filesep 'Turb.mat'])

	 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
	         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
	 
            [ax, ~] = create_axes(fig, 2, 1, 0);
            col = get(groot,'DefaultAxesColorOrder');
            
            xl = Peps.time([1 end]);
            
            %Navg = round(length(Peps.time)/300);
            Navg = round(180/3600/24/diff(Peps.time([1 2])));
            if Navg<1
               Navg=1;
            end
            
            a=1;
               plot(ax(a), Peps.time, Peps.eps_nomask, 'color', [.7 .7 .7 1], 'Linewidth', 1);
               plot(ax(a), Peps.time, Peps.eps, 'color', [col(1,:) .5], 'Linewidth', 1);
               plot(ax(a), Turb.(whichTurb).time, Turb.(whichTurb).eps, 'color', [col(2,:) .5], 'Linewidth', 1);
               pj = 1; p(pj) = plot(ax(a), Peps.time, movmean(Peps.eps, Navg, 'omitnan').*(movmean(isnan(Peps.eps),Navg)<.9), 'color', [col(pj,:)*.5 1], 'Linewidth', 2);
               pj = 2; p(pj) = plot(ax(a),  Turb.(whichTurb).time, movmean( Turb.(whichTurb).eps, Navg, 'omitnan'), 'color', [col(2,:)*.7 1], 'Linewidth', 2);
               pj = 3; p(pj) = plot(ax(a),  Turb.([whichTurb '_ic']).time, movmean( Turb.([whichTurb '_ic']).eps, Navg, 'omitnan'), 'color', [0 1 0 1], 'Linewidth', 2);
               

               legend(p, '\epsilon_{pitot}', ['\epsilon_{\chi} ' whichTurb] , ['\epsilon_{\chi} ' whichTurb '_ic'] );
               yl = 10.^[-9 -2];
               ylim(ax(a), yl);
               xlim(ax(a), xl);
               t = text_corner(ax(a), ['\epsilon [m2/s3]'], 1);
               t = text_corner(ax(a), [unit], -2);
               
               
               set(ax(a), 'Yscale', 'log');

            a=2;
               plot(ax(a), Peps.time, Peps.spd, 'color', [col(1,:) .5], 'Linewidth', 1);
               plot(ax(a), Peps.time,  movmean(Peps.spd, Navg, 'omitnan'), 'color', [col(1,:)*.5 1], 'Linewidth', 2);
               t = text_corner(ax(a), ['speed [m/s]'], 1);

               xlim(ax(a), xl);
               datetick(ax(a), 'keeplimits');
               
                     
               linkaxes(ax, 'x');



          fig1 = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
                  'Papersize',[30 10],'PaperPosition',[0 0 30 10])

                 for a=1:3
                    ax(a) = subplot( 1, 3, a);
                    hold(ax(a), 'on');
                 end

                 a =1;
                 X   = log10(Peps.eps);
                 Y   = log10(clever_interp( Turb.(whichTurb).time, Turb.(whichTurb).eps, Peps.time));
                    XL = '\epsilon_{pitot}';
                    YL = ['\epsilon_{\chi} ' whichTurb];
                    xl = [-10 -3];
                    yl = [min(Y) max(Y)];
                    yl =  xl;
                    binx = [xl(1):(diff(xl)/50):xl(2)];
                    biny = [yl(1):(diff(yl)/50):yl(2)];
                 [hist,~,~,~] = hist2d(binx, biny, X, 0, Y, 0, 3);
                 pcolor(ax(a), binx, biny, hist);    shading(ax(a),'flat');
                 %contourf(ax(a),binx,biny,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
                    plot(ax(a), binx, biny,'k', 'Linewidth', 1);
                    load cmap;
                    colormap(ax(a), cmap.chi);
                    xlabel(ax(a), XL)
                    ylabel(ax(a), YL)
                    ylim(ax(a), yl);
                    xlim(ax(a), xl);
                    % corelation
                       if size(X,2)>size(X,1)
                          X = X';
                       end
                       if size(Y,2)>size(Y,1)
                          Y = Y';
                       end
                       [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
                       t = text_corner(ax(a),  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                                ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);
                 
            
                 a =2;
                 X   = log10(Peps.eps);
                 Y   = log10(clever_interp( Turb.([whichTurb '_ic']).time, Turb.([whichTurb '_ic']).eps, Peps.time));
                    XL = '\epsilon_{pitot}';
                    YL = ['\epsilon_{\chi IC} ' whichTurb];
                    xl = [-10 -3];
                    yl = [min(Y) max(Y)];
                    yl =  xl;
                    binx = [xl(1):(diff(xl)/50):xl(2)];
                    biny = [yl(1):(diff(yl)/50):yl(2)];
                 [hist,~,~,~] = hist2d(binx, biny, X, 0, Y, 0, 3);
                 pcolor(ax(a), binx, biny, hist);    shading(ax(a),'flat');
                 %contourf(ax(a),binx,biny,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
                    plot(ax(a), binx, biny,'k', 'Linewidth', 1);
                    load cmap;
                    colormap(ax(a), cmap.chi);
                    xlabel(ax(a), XL)
                    ylabel(ax(a), YL)
                    ylim(ax(a), yl);
                    xlim(ax(a), xl);
                    % corelation
                       if size(X,2)>size(X,1)
                          X = X';
                       end
                       if size(Y,2)>size(Y,1)
                          Y = Y';
                       end
                       [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
                       t = text_corner(ax(a),  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                                ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);

                 a =3;
                 X   = log10(clever_interp( Turb.([whichTurb]).time, Turb.([whichTurb]).eps, Turb.([whichTurb '_ic']).time));
                 Y   = log10(Turb.([whichTurb '_ic']).eps);
                    XL = ['\epsilon_{\chi} ' whichTurb];
                    YL = ['\epsilon_{\chi IC} ' whichTurb];
                    xl = [-10 -3];
                    yl = [min(Y) max(Y)];
                    yl =  xl;
                    binx = [xl(1):(diff(xl)/50):xl(2)];
                    biny = [yl(1):(diff(yl)/50):yl(2)];
                 [hist,~,~,~] = hist2d(binx, biny, X, 0, Y, 0, 3);
                 pcolor(ax(a), binx, biny, hist);    shading(ax(a),'flat');
                 %contourf(ax(a),binx,biny,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
                    plot(ax(a), binx, biny,'k', 'Linewidth', 1);
                    load cmap;
                    colormap(ax(a), cmap.chi);
                    xlabel(ax(a), XL)
                    ylabel(ax(a), YL)
                    ylim(ax(a), yl);
                    xlim(ax(a), xl);
                    % corelation
                       if size(X,2)>size(X,1)
                          X = X';
                       end
                       if size(Y,2)>size(Y,1)
                          Y = Y';
                       end
                       [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
                       t = text_corner(ax(a),  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                                ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);
