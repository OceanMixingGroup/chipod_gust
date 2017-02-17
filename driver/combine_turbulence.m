%     this script is meant to combine all processed turbulence information 
%
%    created:
%     Johannes Becherer
%     Tue Aug 23 12:26:41 PDT 2016

clear all;
close all;

%_____________________set flags______________________
   do_combine  =  1; % do actually the averaging calculation (can take a couple o minutes)
   do_plot     =  1; % generate a comparison plot between the different estimates 

   % if you want to restrict the time range that should be combined use the following
   time_range(1)  = datenum(2000, 1, 1, 0, 0, 0); 
   time_range(2)  = datenum(2030, 1, 1, 0, 0, 0); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________find all available chi data______________________
if(do_combine)
   d = dir([basedir '/proc/chi/']);
   % eliminate directories
   d = d(~[d(:).isdir]);

   % run through all chi files
   for i = 1:length(d)

      mat_test = strfind(d(i).name, '.mat');
      chi_test = strfind(d(i).name, 'chi');
      ic_test  = strfind(d(i).name, 'ic.');

      if ~isempty(mat_test) & ~isempty(chi_test)

         ID = d(i).name(1:mat_test-1);
         disp(['adding ' ID ]);
         load([basedir '/proc/chi/' ID '.mat'])
         
         % find desired time range
         iiTrange = find( chi.time>= time_range(1) & chi.time<= time_range(2) );

         % get list of all fields to average
         ff = fields(chi);

         if isempty(ic_test) % not inertial convective estimate
            %% average data
            ww =  round(600*diff(chi.time(1:2))*3600*24); % averaging window

            for f = 1:length(ff)  % run through all fields in chi
               if ( length(chi.(ff{f})) == length(chi.time) )

                  % deglitch chi and eps
                  if strcmp(ff{f},'eps') | strcmp(ff{f},'chi')
                     chi.(ff{f}) = deglitch(chi.(ff{f}), ww, 2,'b');
                  end


                  Turb.(ID).(ff{f}) = moving_average( chi.(ff{f})(iiTrange), ww, ww );
               end
            end

         else % IC estimate is already time averaged

            for f = 1:length(ff)  % run through all fields in chi
               if ( length(chi.(ff{f})) == length(chi.time) )
                  Turb.(ID).(ff{f}) = chi.(ff{f})(iiTrange);
               end
            end

         end

      end
   end

%_____________________save combined structure______________________
   save([basedir '/proc/Turb.mat'], 'Turb');

    
end

%_____________________comparison plot______________________
if do_plot
   
   load([basedir '/proc/Turb.mat']);
   ff = fields(Turb);

    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[30 20],'PaperPosition',[0 0 30 20])
    
         [ax, ~] = create_axes(fig, 4, 1, 0);
      
         col = get(groot,'DefaultAxesColorOrder');
         
         a=1;
         for f = 1:length(ff)
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).chi, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['\chi [K^2/s]'], 1);
         set(ax(a), 'Yscale', 'log');
         yl = [1e-10 1e-2];
         ylim(ax(a), yl);
         
         a=2;
         for f = 1:length(ff)
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).eps, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['\epsilon [m^2/s^3]'], 1);
         set(ax(a), 'Yscale', 'log');
         yl = [1e-10 1e-2];
         ylim(ax(a), yl);
         
         
         a=3;
         for f = 1:length(ff)
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).N2, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['N^2 [s^{-2}]'], 1);
         set(ax(a), 'Yscale', 'log');
         yl = [1e-6 1e-2];
         ylim(ax(a), yl);
         

         a=4;
         for f = 1:length(ff)
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).spd, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['|u| [m/s]'], 1);
         datetick(ax(a), 'keeplimits');
         yl = [0 2];
         ylim(ax(a), yl);

         linkaxes(ax, 'x');
         
            
   %---------------------histogram plots----------------------
      squeeze_axes(ax, .8, 1)

      pos = get(ax(1), 'Position');
      axh(1) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.01 pos(4)] ) 
      
      pos = get(ax(2), 'Position');
      axh(2) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.01 pos(4)] ) 

         a=1;
         hold(axh(a),'all');
         yl = log10(get(ax(a), 'Ylim'));
         bins = yl(1):diff(yl)/100:yl(2);
         for f = 1:length(ff)
            [Nchi,~] = histcounts( log10(Turb.(ff{f}).chi) , bins);
            pj = f; p(pj) = plot(axh(a), Nchi , bins(1:end-1)+diff(bins(1:2)*.5), 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         ylim(axh(a), yl);
         set(axh, 'Yticklabel', {}, 'Xticklabel', {})
      
         a=2;
         hold(axh(a),'all');
         yl = log10(get(ax(a), 'Ylim'));
         bins = yl(1):diff(yl)/100:yl(2);
         for f = 1:length(ff)
            [Nchi,~] = histcounts( log10(Turb.(ff{f}).eps) , bins);
            pj = f; p(pj) = plot(axh(a), Nchi , bins(1:end-1)+diff(bins(1:2)*.5), 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         ylim(axh(a), yl);
         set(axh, 'Yticklabel', {}, 'Xticklabel', {})
         

   %---------------------legend----------------------
      pos = get(ax(3), 'Position');
      axl = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.01 pos(4)] ) 
      hold(axl,'on');
      

         for f = 1:length(ff)
            pj = f; p(pj) = plot(axl, [0 1] ,[0 1], 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         legend(p, ff);
         set(axl, 'visible', 'off')
         ylim(axl,[-1 -.5])


   %---------------------save imagage----------------------

   print(gcf,'../pics/Compare_Turb.png','-dpng','-r200','-painters')
   
end