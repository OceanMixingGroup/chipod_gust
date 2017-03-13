% Reads and processes PMELs crazy time vectors into MATLAB format
%           [timeout] = ProcessPmelTime(fname)
% fname is valid netCDF file name with path etc.
% returns time vector in MATLAB format
function [timeout] = ProcessPmelTime(fname)

    % thanks johannes
    timein = double(ncread(fname, 'time'));
    timeFormat = ncreadatt(fname, 'time', 'units');

    switch timeFormat(1:3)
      case 'day'
        timeUnit = 1;
      case 'hou'
        timeUnit = 1/(24);
      case 'min'
        timeUnit = 1/(60*24);
      case 'sec'
        timeUnit = 1/(3600*24);
      otherwise
        disp('unknown Time unit.....check NETCDF file')
    end
    % time off set
    timeOffset = datenum(timeFormat(end-18:end));

    timeout = timein*timeUnit + timeOffset ;
end
