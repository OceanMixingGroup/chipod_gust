% Return temp, salinity file names
function fname = GetPmelFileName(varname, lon, lat, datadir, ...
                                 array, res)

    if lon > 0
        lonstr = 'e';
    else
        lonstr = 'w';
    end

    if lat < 0
        latstr = 's';
    else
        latstr = 'n';
    end

    if ~strcmpi(res, 'dy')
        subdir = ['high_resolution/' res '/'];
        fname = [datadir subdir varname   num2str(abs(lat)), ...
                 latstr,num2str(abs(lon)),lonstr,'_' res '.cdf'];
    else
        subdir = [upper(array) '/'];
        fname = [datadir subdir varname '_xyzt_dy.cdf'];
    end

    fname = ls(fname);
    fname = fname(1:end-1); % thoo. ls tacks on an extra space at the end
end
