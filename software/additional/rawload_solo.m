function [T] = rawload_solo(fid)
%%  [T] = rawload_solo(fid)
%
%  This function reads in the raw files (sqlite, rsk)
%  from SoloTs and returns a  mat struction structure 
%
%  INPUT
%     fid      :  path to raw-file
%
%  OUTPUT
%     T.time   :  matlab time
%     T.T      :  temperature
%
%
%  !!!! NOTE this function requirres MATLAB2016a or newer!!!!
%
%   created by: 
%        Johannes Becherer
%        Sat Jul  8 00:13:24 GMT 2017


if nargin<1
    [raw_name,rawdir]=uigetfile('*.*','Load Binary File');
    fid=[rawdir raw_name];
    if raw_name==0
        error('File not found')
        return
    end
end


% open data base file
db1 = sqlite(fid);

% get data block
data = fetch(db1,'select * from data');

% convert unix time from t logger into matlab time
timetmp = datenum(1970,1,1,0,0,double(cell2mat(data(:,1))/1000));
[T.time, ii]  = sort(timetmp); % the data require chronological sorting

% temperature
ttmp    = double(cell2mat(data(:,2)));
T.T     = ttmp(ii);

