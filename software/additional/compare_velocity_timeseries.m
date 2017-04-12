function [fig] = compare_velocity_timeseries(a_time, a_U, a_L, p_time, p_U, p_L, vis)
%%    [fig] = compare_velocity_timeseries(a_time, a_U , a_L, p_time, p_U, p_L, vis)
%           
%           Generates a monster plot that compares both velocity time series
%
%           INPUT
%              a_time      :  (adcp time time vector)
%              a_U         :  complex velocity for a
%              a_L         :  Label for a
%              r_time      :  (pitot time time vector)
%              r_U         :  complex velocity for p
%              p_L         :  Label for p
%              vis         :  figure visibility (default 'on')
%
%   created by: 
%        Johannes Becherer
%        Thu Nov 10 10:26:01 PST 2016


if nargin<7
   vis = 'on';
end

%_____________________identify common limits______________________
xl = [max([a_time(1), p_time(1)]) min([a_time(end), p_time(end)])];
aii = find( a_time>=xl(1) & a_time<xl(end) );
pii = find( p_time>=xl(1) & p_time<xl(end) );

%_____________________common time vector______________________
% choose the courser time resolution
adt = diff(a_time(1:2));
pdt = diff(p_time(1:2));
if adt>pdt
   C.time = a_time(aii);
   C.a_U  = a_U(aii);
   C.p_U  = interp1( p_time(pii), qbutter(p_U(pii), pdt/adt), C.time);
   
else
   C.time = p_time(pii);
   C.p_U  = p_U(pii);
   C.a_U  = interp1( a_time(aii), qbutter(a_U(aii), adt/pdt), C.time);
end

%_____________________Plot figure______________________
 col = get(groot,'DefaultAxesColorOrder');

 % limits
 sl = [0 max([ max(abs(a_U)), max(abs(p_U)) ])]; 

 fig = figure('Color',[1 1 1],'visible',vis,'Paperunits','centimeters',...
         'Papersize',[15 20],'PaperPosition',[0 0 15 20])
   
         % angle axes
          [axa, ~] = create_axes(fig,1,2, 0);
            squeeze_axes(axa, 1,.3);
            shift_axes(axa, 0.01,-.03);
          
         % time sereis axes
          [axt, ~] = create_axes(fig,1,1, 0);
            squeeze_axes(axt, 1,.25);
            shift_axes(axt, 0.01,.64);
          
         % histo axes
          [axh, ~] = create_axes(fig,1,2, 0);
            squeeze_axes(axh, 1,.3);
            shift_axes(axh, 0.01,.32);

   %---------------------plot time series----------------------
   ax = axt; a=1;
   hold(ax(a),'on');
   plot(ax(a), xl, [0 0], 'k--', 'Linewidth', 1);
   ww = ceil(length(C.time)/100);
   pj = 1; p(pj) = plot(ax(a), C.time, real(C.a_U), 'color', [col(pj,:) .3], 'Linewidth', 1);
   pj = 2; p(pj) = plot(ax(a), C.time, imag(C.a_U), 'color', [[0 0 1]./1.5 .3], 'Linewidth', 1);
   pj = 3; p(pj) = plot(ax(a), C.time, real(C.p_U), 'color', [col(2,:)  .3], 'Linewidth', 1);
   pj = 4; p(pj) = plot(ax(a), C.time, imag(C.p_U), 'color', [[1 0 0]./1.5 .3], 'Linewidth', 1);

   pj = 1; p(pj) = plot(ax(a), moving_average(C.time, ww, ww), moving_average(real(C.a_U), ww, ww),'color', [col(pj,:) .9], 'Linewidth', 2);
   pj = 2; p(pj) = plot(ax(a), moving_average(C.time, ww, ww), moving_average(imag(C.a_U), ww, ww),'color', [[0 0 1]./1.5 .9], 'Linewidth', 2);
   pj = 3; p(pj) = plot(ax(a), moving_average(C.time, ww, ww), moving_average(real(C.p_U), ww, ww),'color', [col(2,:)  .9], 'Linewidth', 2);
   pj = 4; p(pj) = plot(ax(a), moving_average(C.time, ww, ww), moving_average(imag(C.p_U), ww, ww),'color', [[1 0 0]./1.5 .9], 'Linewidth', 2);
      xlim(ax(a), xl);

      ylim(ax(a), [-sl(2) sl(2)]);
      
      set(ax(a), 'Xtick', [ceil(xl(1)):ceil(diff(xl)/5):floor(xl(2))])
      datetick(ax(a), 'keeplimits', 'keepticks');
      %xlabel(ax(a), datestr(nanmedian(C.time), 'dd mmm yyyy'));
      legend(p, ['u_{' a_L '}'], ['v_{' a_L '}'], ['u_{' p_L '}'], ['v_{' p_L '}'],...
            'orientation', 'horizontal', 'location', 'southwest');
      ylabel(ax(a), 'm s^{-1}');

   %---------------------angle plot----------------------
   a=2
   hold(axa(a),'off')
   rose(axa(a), angle(C.a_U), 50);
   hold(axa(a),'on')
   rose(axa(a), angle(C.p_U), 50);
      t = text_corner(axa(a), [a_L], -7);
      t.Color = col(1,:);
      t = text_corner(axa(a), {'';p_L},-7);
      t.Color = col(2,:);
  t = text_corner(axa(2), ['flow direction'], -9);
   
   
   %---------------------angle histogram----------------------
  %X{1}   = angle(C.a_U);
  %X{2}   = angle(C.p_U);
  %al = [-pi pi];
  %bins = [al(1):(diff(al)/50):al(2)];
  %for i =1:length(X);
  %   [Ncnt{i},~] = histcounts( X{i} , bins);
  %   pj = i; p(pj) = plot( axa(2), bins(1:end-1)+diff(bins(1:2)*.5), ...
  %                        Ncnt{i} , 'color', [col(pj,:) 1], 'Linewidth', 2);   
  %   patch( [bins(1) bins(1:end)]+diff(bins(1:2)*.5), [0 Ncnt{i} 0], ...
  %            col(i,:)*.5+.5,'Facealpha',.3, 'Parent', axa(2));
  %end
  %xlabel(axa(2), ['[rad]'])
  %   t = text_corner(axa(2), [a_L], 3);
  %   t.Color = col(1,:);
  %   t = text_corner(axa(2), {'';p_L}, 3);
  %   t.Color = col(2,:);
  %set(axa(2), 'Ycolor', [1 1 1]);
  %xlim(axa(2), al);
  %t = text_corner(axa(2), ['flow direction'], 1);
  
  %_____________________speed ratio______________________
  a=1;
  X{1}   = log10(abs(C.p_U( abs(C.a_U)>.05 &  abs(C.p_U)>.05))./abs(C.a_U(abs(C.a_U)>.05 &  abs(C.p_U)>.05)));
  sl = [min(X{1}) max(X{1})];
  bins = [sl(1):(diff(sl)/50):sl(2)];
  for i =1;
     [Ncnt{i},~] = histcounts( X{i} , bins);
     pj = i; p(pj) = plot( axa(a), bins(1:end-1)+diff(bins(1:2)*.5), ...
                          Ncnt{i} , 'color', [col(4,:) 1], 'Linewidth', 2);   
     patch( [bins(1) bins(1:end)]+diff(bins(1:2)*.5), [0 Ncnt{i} 0], ...
              col(4,:)*.5+.5,'Facealpha',.3, 'Parent', axa(a));
  end
  plot(axa(a), [0 0], [0 max(Ncnt{1})], '--k', 'Linewidth', 1);
  
  xlabel(axa(a), ['log_{10}|u_{' p_L '}|/|u_{' a_L '}|'])
  t = text_corner(axa(a), ['avg ' num2str(nanmean(X{1}), '%1.2f') ], 3);
  t.Color = col(4,:);
  t = text_corner(axa(a), {'';['med ' num2str(nanmedian(X{1}), '%1.2f') ]}, 3);
  t.Color = col(4,:);
  set(axa(a), 'Ycolor', [1 1 1]);
  xlim(axa(a), sl);
  
   

   
   %---------------------speed histogram----------------------
   X{1}   = abs(C.a_U);
   X{2}   = abs(C.p_U);
   sl = [min(X{1}) max(X{1})];
   bins = [sl(1):(diff(sl)/50):sl(2)];
   for i =1:length(X);
      [Ncnt{i},~] = histcounts( X{i} , bins);
      pj = i; p(pj) = plot( axh(2), bins(1:end-1)+diff(bins(1:2)*.5), ...
                           Ncnt{i} , 'color', [col(pj,:) 1], 'Linewidth', 2);   
      patch( [bins(1) bins(1:end)]+diff(bins(1:2)*.5), [0 Ncnt{i} 0], ...
               col(i,:)*.5+.5,'Facealpha',.3, 'Parent', axh(2));
   end
   xlabel(axh(2), ['[m s^{-1}]'])
   set(axh(2), 'Ycolor', [1 1 1]);
   xlim(axh(2), sl);
   t = text_corner(axh(2), ['speed'], 1);
      t = text_corner(axh(2), ['    avg   med '], 3);
      t = text_corner(axh(2), {'' ;[a_L '  ' num2str(nanmean(X{1}), '%1.2f') '  ' num2str(nanmedian(X{1}), '%1.2f') ' ']}, 3);
      t.Color = col(1,:);
      t = text_corner(axh(2), {'' ;'' ;[p_L '  ' num2str(nanmean(X{2}), '%1.2f') '  ' num2str(nanmedian(X{2}), '%1.2f')]}, 3);
      t.Color = col(2,:);

   
   %_____________________2D hist______________________
   bins = [0:(diff(sl)/50):sl(2)];
   [hist,mn,mdn,md] = hist2d(bins, bins, abs(C.a_U), 0, abs(C.p_U), 0, 3);
   %pcolor(axh(1), bins, bins, hist);
   %    shading(axh(1),'flat');
   contourf(axh(1),bins,bins,hist,[0:.1:1]*max(max(hist)), 'edgecolor', 'none');
       load cmap;
      colormap(axh(1), cmap.chi);
       colormap(axh(1), cmap.chi);
       plot(axh(1), bins, bins,'k', 'Linewidth', 1);
      ylabel(axh(1), ['|u_{' p_L '}| [m s^{-1}]'])
      xlabel(axh(1), ['|u_{' a_L '}| [m s^{-1}]'])
            % corelation
            X   = abs(C.a_U);
            Y   = abs(C.p_U);
               if size(X,2)>size(X,1)
                  X = X';
               end
               if size(Y,2)>size(Y,1)
                  Y = Y';
               end
               [r, ~, rL, rH] = corrcoef(X( ~isnan(X) & ~isnan(Y) ), Y( ~isnan(X) & ~isnan(Y) ));
               t = text_corner(axh(a),  {['r = ' num2str(r(2)*100, '%2.1f')  ' % ']; ...
                           ['  [ ' num2str(rL(2)*100, '%2.1f') ', ' num2str(rH(2)*100, '%2.1f') ' ]']}, 6);
               
            
       
   
   %_____________________abc______________________
   abc='abcdefghijklmnopqrst';
   ax(1) = axt(1);
   ax(2) = axh(1);
   ax(3) = axh(2);
   ax(4) = axa(1);
   ax(5) = axa(2);
   for a = 1:(size(ax,1)*size(ax,2))
      t = text_corner(ax(a), abc(a), 9);
      t.BackgroundColor = [1 1 1 .7];
      set(ax(a), 'Tickdir', 'out')
   end
      set(axt, 'box', 'on')
      set(axh(1), 'box', 'on')
   
   
   
      
   
   
   
   
   
   



