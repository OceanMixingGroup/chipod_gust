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



   % correlation

   % bring on common time step

   Peps.dt = diff(Peps.time([1 2]));
   Dt = 3*60/24/3600;
   NP60sec = round((Dt)/Peps.dt);
   C.time         = Peps.time([1:NP60sec:end]);
   C.peps         = interp1( Peps.time, movmean(Peps.eps, NP60sec, 'omitnan'), C.time);
   C.eps_bbl      = interp1( Peps.time, abs((2.2e-3)^1.5/(.4*5)*movmean(Peps.spd, NP60sec, 'omitnan').^3), C.time);
   C.eps_chi      = interp1( Turb.(whichTurb).time, movmean(Turb.(whichTurb).eps, round(Dt/diff(Turb.(whichTurb).time([1 2]))), 'omitnan'), C.time); 
   C.eps_chi_ic   = interp1( Turb.([whichTurb '_ic']).time, movmean(Turb.([whichTurb '_ic']).eps, ceil(Dt/diff(Turb.([whichTurb '_ic']).time([1 2]))), 'omitnan'), C.time); 

   cor_ww  = 30; % 3xmin
   C.cor_peps_vc  = movcorr( C.peps', C.eps_chi', cor_ww, 'omitnan' );
   C.cor_peps_ic  = movcorr( C.peps', C.eps_chi_ic', cor_ww, 'omitnan' );
   C.cor_vc_ic    = movcorr( C.eps_chi', C.eps_chi_ic', cor_ww, 'omitnan' );
   C.cor_peps_bbl = movcorr( C.eps_chi', C.eps_bbl', cor_ww, 'omitnan' );

    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[30 20],'PaperPosition',[0 0 30 20])

            col = get(groot,'DefaultAxesColorOrder');
    
            [ax, ~] = create_axes(fig, 4, 1, 0);
            
            a=1;
            plot(ax(a), C.time, C.peps, 'Linewidth', 1);
            plot(ax(a), C.time, C.eps_chi, 'Linewidth', 1);
            plot(ax(a), C.time, C.eps_chi_ic, 'Linewidth', 1);
            plot(ax(a), C.time, C.eps_bbl, 'Linewidth', 1);
            set(ax(a), 'Yscale', 'log');

            a=2;
            plot(ax(a), C.time, C.cor_peps_vc, 'Linewidth', 1);
            plot(ax(a), C.time, C.cor_peps_ic, 'Linewidth', 1);
            plot(ax(a), C.time, C.cor_vc_ic, 'Linewidth', 1);
            plot(ax(a), C.time, C.cor_peps_bbl, 'Linewidth', 1);

            
            a=3;
               plot(ax(a), Peps.time, Peps.spd, 'color', [col(1,:) .5], 'Linewidth', 1);
               plot(ax(a), Peps.time,  movmean(Peps.spd, NP60sec, 'omitnan'), 'color', [col(1,:)*.5 1], 'Linewidth', 2);
               t = text_corner(ax(a), ['speed [m/s]'], 1);
               
            a=4;
               plot(ax(a),Turb.(whichTurb).time, movmean( Turb.(whichTurb).N2, 4, 'omitnan'), 'color', [col(1,:) .5], 'Linewidth', 1);
               t = text_corner(ax(a), ['N2 [1/s2]'], 1);
               set(ax(a), 'Yscale', 'log');

               datetick(ax(a), 'keeplimits');
               linkaxes(ax, 'x');
               
                  


	 fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
	         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
	 
            [ax, ~] = create_axes(fig, 3, 1, 0);
            col = get(groot,'DefaultAxesColorOrder');
            
            xl = Peps.time([1 end]);
            
            %Navg = round(length(Peps.time)/300);
            Navg = round(120/3600/24/diff(Peps.time([1 2])));
            if Navg<1
               Navg=1;
            end
            
            a=1;
               plot(ax(a), Peps.time, Peps.eps_nomask, 'color', [.7 .7 .7 1], 'Linewidth', 1);
               plot(ax(a), Peps.time, Peps.eps, 'color', [col(1,:) .5], 'Linewidth', 1);
               plot(ax(a), Turb.(whichTurb).time, Turb.(whichTurb).eps, 'color', [col(2,:) .5], 'Linewidth', 1);
               pj = 1; p(pj) = plot(ax(a), Peps.time, movmean(Peps.eps, Navg, 'omitnan').*(movmean(isnan(Peps.eps),Navg)<.9), 'color', [col(pj,:)*.5 1], 'Linewidth', 2);
               pj = 2; p(pj) = plot(ax(a),  Turb.(whichTurb).time, movmean( Turb.(whichTurb).eps, 4, 'omitnan'), 'color', [col(2,:)*.7 1], 'Linewidth', 2);
               if isfield(Turb,[whichTurb '_ic'] )
                  pj = 3; p(pj) = plot(ax(a),  Turb.([whichTurb '_ic']).time, movmean( Turb.([whichTurb '_ic']).eps, 1, 'omitnan'), 'color', [0 1 0 1], 'Linewidth', 2);
                  legend(p, '\epsilon_{pitot}', ['\epsilon_{\chi} ' whichTurb] , ['\epsilon_{\chi} ' whichTurb '_ic'] );
               else
                  legend(p, '\epsilon_{pitot}', ['\epsilon_{\chi} ' whichTurb] );
               end
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
               
            a=3;
               plot(ax(a),Turb.(whichTurb).time, movmean( Turb.(whichTurb).N2, 4, 'omitnan'), 'color', [col(1,:) .5], 'Linewidth', 1);
               t = text_corner(ax(a), ['N2 [1/s2]'], 1);
               set(ax(a), 'Yscale', 'log');

               xlim(ax(a), xl);
               datetick(ax(a), 'keeplimits');
                     
               linkaxes(ax, 'x');



          fig1 = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
                  'Papersize',[30 10],'PaperPosition',[0 0 30 10])

                 for a=1:4
                    ax(a) = subplot( 2, 2, a);
                    hold(ax(a), 'on');
                 end

                 a =1;
                 X   = log10(C.peps);
                 Y   = log10(C.eps_chi);
                 %Y(C.cor_peps_vc<.5 ) = nan;
                 %Y(Y<-7.5) = nan;
                 %X = X(~isnan(Y));
                 %Y = Y(~isnan(Y));
                 
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
                    plot(ax(a), binx, biny-1,'k', 'Linewidth', 1);
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
                 if isfield(Turb,[whichTurb '_ic'] )
                 X   = log10(C.peps);
                 Y   = log10(C.eps_chi_ic);
                    XL = '\epsilon_{pitot}';
                    YL = ['\epsilon_{\chi IC} ' whichTurb];
                    xl = [-10 -3];
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
                 X   = log10(C.eps_chi);
                 Y   = log10(C.eps_chi_ic);
                    XL = ['\epsilon_{\chi} ' whichTurb];
                    YL = ['\epsilon_{\chi IC} ' whichTurb];
                    xl = [-10 -3];
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


                a=4; 
                 X   = log10(C.peps);
                 Y   = log10(abs(C.eps_bbl));
                    XL = ['\epsilon_{pitot} ' whichTurb];
                    YL = ['\epsilon_{bbl} ' whichTurb];
                    xl = [-10 -3];
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
                 end
