function [out] = choose_color(name,outstr)
% [out] = choose_color(name,outstr)
%
% When processing chipods, there can be a lot of different lines. Here,
% we're tring to define the colors for each processing mode, so that they
% will always be plotted with the same colors. I have chosen 16
% distinguishable colors:
%
% mi11 - red            mi11_ic - pink
% mi22 - purple         mi22_ic - lilac
% mm1 - blue            mm1_ic - light blue
% mm2 - dark green      mm2_ic - light green
% pi11 - black          pi11_ic - gray
% pi22 - brown          pi22_ic - tan
% pm1 - orange          pm1_ic - coral
% pm2 - dark yellow     pm2_ic - light yellow
%
% name   = name of processing mode (e.g. mi11, pm1, etc). must be a string
%
% outstr = 'color'    to return 3 digit rgb vector
% outstr = 'style'    to return linestyle ('-', ':', '--')
% outstr = 'width'    to return linewidth


% define possible colors
cols =      [0         0         0;...          % black             pi11
             0.9333    0.1098    0.1412;...     % red               mi11
             0         0.5000    1.0000;...     % blue              mm1   
             0.6667    0.6667    0.6667;...     % gray              pi11_ic
             0.5725    0.1529    0.5608;...     % purple            mi22
             0.2330    0.5500    0.0900;...     % dark green        mm2
             0.9725    0.5804    0.1137;...     % orange            pm1
             0.7804    0.7333    0.0039;...     % visible yellow    pm2
             0.6510    0.4863    0.3216;...     % brown             pi22
             0.7451    0.6392    0.8549;...     % lilac             mi22_ic
             0.9608    0.6745    0.8000;...     % pink              mi11_ic
             0.4660    0.8000    0.1880;...     % light green       mm2_ic  
             0.6000    0.8000    0.9800;...     % light blue        mm1_ic
             0.8255    0.7431    0.6608;...     % tan               pi22_ic
             0.9863    0.7902    0.5568;...     % coral             pm1_ic  
             0.8902    0.8666    0.5020];       % light yellow      pm2_ic


% define styles
if strcmp(name(end-1:end),'ic')
    style = '-';
    width = 1;
else
    style = '-';
    width = 1;
end

% find the color that relates to each input
if strcmp(name(5:end),'pi11')
    color = cols(1,:);
elseif strcmp(name(5:end),'mi11')
    color = cols(2,:);
elseif strcmp(name(5:end),'mm1')
    color = cols(3,:);
elseif strcmp(name(5:end),'pi11_ic')
    color = cols(4,:);
elseif strcmp(name(5:end),'mi22')
    color = cols(5,:);
elseif strcmp(name(5:end),'mm2')
    color = cols(6,:);
elseif strcmp(name(5:end),'pm1')
    color = cols(7,:);
elseif strcmp(name(5:end),'pm2')
    color = cols(8,:);
elseif strcmp(name(5:end),'pi22')
    color = cols(9,:);
elseif strcmp(name(5:end),'mi22_ic')
    color = cols(10,:);
elseif strcmp(name(5:end),'mi11_ic')
    color = cols(11,:);
elseif strcmp(name(5:end),'mm2_ic')
    color = cols(12,:);
elseif strcmp(name(5:end),'mm1_ic')
    color = cols(13,:);
elseif strcmp(name(5:end),'pi22_ic')
    color = cols(14,:);
elseif strcmp(name(5:end),'pm1_ic')
    color = cols(15,:);
elseif strcmp(name(5:end),'pm2_ic')
    color = cols(16,:);
else 
    disp('repeating colors for mi112, mi212, pi112 or pi212')
    if strcmp(name(5:end),'mi112')
        color = cols(2,:);
    elseif strcmp(name(5:end),'mi112_ic')
        color = cols(11,:);
    elseif strcmp(name(5:end),'mi212')
        color = cols(5,:);
    elseif strcmp(name(5:end),'mi212_ic')
        color = cols(10,:);
    elseif strcmp(name(5:end),'pi112')
        color = cols(1,:);
    elseif strcmp(name(5:end),'pi112_ic')
        color = cols(4,:);
    elseif strcmp(name(5:end),'pi212')
        color = cols(9,:);
    elseif strcmp(name(5:end),'pi212_ic')
        color = cols(15,:);
    end
end

% determine output
if strcmp(outstr,'color')
    out = color;
elseif strcmp(outstr,'style')
    out = style;
elseif strcmp(outstr,'width')
    out = width;
end
    
    
    
end

