% plot_Turb.m
%
% load Turb.m and plot all of the variables that have been processed

%%%%%%%%%%%%%%%%%%% DEFINE PARAMETERS %%%%%%%%%%%%%%%%%%%

% SAVE FIG
% The fig file can be very large and take a huge amount of time to save.
% Only turn on with smaller datasets
do_save_fig = 0;

% RUNNAME
% Name of folder in which Turb.mat is saved (folder within ../proc/)
% For example, if Turb.mat is in ../proc/Turb_min_dTdz_1e-3_60sec/Turb.mat,
% then runname = 'Turb_min_dTdz_1e-3_60sec'
% If Turb.mat is saved in ../proc/, leave runname = NaN;
runname = NaN;

% PROCESSING MODES TO PLOT
% There are up to 24 processing modes, which would you like to plot?
% Note: if any of these have not been processed, they will not plot
% automatically.
if ~exist('do', 'var')
    do = ChooseEstimates();
end

% filter out estimates here:
% EXAMPLE 1: do = ChooseEstimates(do, 'no_ic')
% EXAMPLE 2: do = ChooseEstimates([], 'none'); do.chi_pm1 = 1; % only pm1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set directories 
	addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines
    here    =   pwd;                % mfiles folder
    basedir =   here(1:(end-6));    % substract the mfile folder
    savedir =   [basedir 'proc/'];  % directory directory to save data
    unit    = chi_get_unit_name(basedir); % get unit name

% load in Turb.mat
tic
if ~isnan(runname)
    disp(['loading ../proc/' runname '/Turb.mat'])
	load([basedir '/proc/' runname '/Turb.mat']);
else
    disp('loading ../proc/Turb.mat')
	load([basedir '/proc/Turb.mat']);
end
toc

ff = fields(Turb);
ff = {ff{1:end-1}}'; % remove readme structure

fig = CreateFigure;
[ax, ~] = create_axes(fig, 4, 1, 0);
      
 
%%%%%%%%%%%%%%%%% plot time series %%%%%%%%%%%%%%%%%
for a = 1:4
    % define parameters for each subplot
    if a == 1
        var = 'chi';
        labelstr = '\chi [K^2/s]';
        yl = [1e-10 1e-4];
        yscale = 'log';
    elseif a == 2
        var = 'eps';
        labelstr = '\epsilon [m^2/s^3]';
        yl = [1e-10 1e-4];
        yscale = 'log';
    elseif a == 3
        var = 'dTdz';
        labelstr = 'dTdz [C/m]';
        yl = [-1e-1, 0.3]; 
        yscale = 'linear';
    elseif a == 4
        var = 'spd';
        labelstr = '|u| [m/s]';
        yl = [0 2];
        yscale = 'linear';
    end
    
    % plot
    for f = 1:length(ff)
        if ~isstruct(Turb.(ff{f})), continue; end
        if do.(ff{f}) == 1
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).(var),...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'));
        end
    end
    t = text_corner(ax(a), labelstr, 1);
    set(ax(a), 'Yscale', 'log');    
    ylim(ax(a), yl);
    set(ax(a),'tickdir','out')
    datetick(ax(a), 'keeplimits');
    
    % add some extras
    if a == 1 | a == 2
        set(ax(a),'YTick',10.^(-10:1:0))
    end
    if a == 3
        set(ax(a), 'XTickLabels', [])
        plot(ax(a), xlim ,[0 0], 'k-')
        axes(ax(a))
        symlog('y', -3);
    end
    if a ~= 4
        set(ax(a),'XTickLabel',{' '})
    end
    set(ax(4),'YScale','linear')

end

linkaxes(ax, 'x');


%%%%%%%%%%%%%%%%% histograms %%%%%%%%%%%%%%%%%

ms = 8; %marker size of means and medians

squeeze_axes(ax, .8, 1)

pos = get(ax(1), 'Position');
axh(1) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.04 pos(4)] );
      
pos = get(ax(2), 'Position');
axh(2) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.04 pos(4)] );

a=1;
hold(axh(a),'all');
yl = log10(get(ax(a), 'Ylim'));
bins = yl(1):diff(yl)/100:yl(2);
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1, continue; end
    if do.(ff{f}) == 1
        [Nchi,~] = histcounts( log10(Turb.(ff{f}).chi) , bins);
        pj = f; p(pj) = plot(axh(a), Nchi , bins(1:end-1)+diff(bins(1:2)*.5),...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'));
    end
end
xlims = get(axh(a),'XLim');
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1, continue; end
    if do.(ff{f}) == 1 & isfield(Turb.(ff{f}), 'stats')
        plot(axh(a), xlims(2)*0.9 ,log10( Turb.(ff{f}).stats.chimean),'<',...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'),'markersize',ms);
        plot(axh(a), xlims(2)*0.95 ,log10( Turb.(ff{f}).stats.chimedian),'+',...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'),'markersize',ms);
    end
end
ylim(axh(a), yl);
set(axh(a), 'Xticklabel', {},'tickdir','out','YAxisLocation','right',...
    'YTick',-10:1:0)
xl = get(axh(a),'XLim');
text(axh(a),xl(2)*0.95,yl(2),' median -->','rotation',-90,'fontsize',10,...
    'horizontalalignment','left','verticalalignment','middle');
text(axh(a),xl(2)*0.9,yl(2),' mean -->','rotation',-90,'fontsize',10,...
    'horizontalalignment','left','verticalalignment','middle');
      
a=2;
hold(axh(a),'all');
yl = log10(get(ax(a), 'Ylim'));
bins = yl(1):diff(yl)/100:yl(2);
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1, continue; end
    if do.(ff{f}) == 1
        [Nchi,~] = histcounts( log10(Turb.(ff{f}).eps) , bins);
        pj = f; p(pj) = plot(axh(a), Nchi , bins(1:end-1)+diff(bins(1:2)*.5),...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'));
    end
end
xlims = get(axh(a),'XLim');
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1, continue; end
    if do.(ff{f}) == 1  & isfield(Turb.(ff{f}), 'stats')
        plot(axh(a), xlims(2)*0.9 , log10(Turb.(ff{f}).stats.epsmean),'<',...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'),'markersize',ms);
        plot(axh(a), xlims(2)*0.95 , log10(Turb.(ff{f}).stats.epsmedian),'+',...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'),'markersize',ms);
    end
end
ylim(axh(a), yl);
% set(axh, 'Yticklabel', {}, 'Xticklabel', {},'tickdir','out')
set(axh(a), 'Xticklabel', {},'tickdir','out','YAxisLocation','right',...
    'YTick',-10:1:0)



%%%%%%%%%%%%%%%%% legend %%%%%%%%%%%%%%%%%
pos = get(ax(3), 'Position');
axl = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.04 pos(4)] );
hold(axl,'on');

clear nleg p ffstr
nleg = 1;
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1, continue; end
    if do.(ff{f}) == 1
        pj = nleg; p(pj) = plot(axl, [0 1] ,[0 1],...
                'color', choose_color(ff{f}(5:end),'color'), ...
                'LineWidth', choose_color(ff{f}(5:end),'width'), ...
                'LineStyle', choose_color(ff{f}(5:end),'style'));
        ffstr{nleg} = fix_underscore(ff{f}(5:end));
        nleg = nleg+1;
    end
end
legend(p, ffstr);
set(axl, 'visible', 'off')
ylim(axl,[-1 -.5])


%%%%%%%%%%%%%%%%% save %%%%%%%%%%%%%%%%%

print(gcf,'../pics/Compare_Turb.png','-dpng','-r200','-painters')
if do_save_fig == 1
    disp('saving ../pics/Compare_Turb.fig')
    tic
    savefig(gcf,'../pics/Compare_Turb.fig')
    toc
end




   
