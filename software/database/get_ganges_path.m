% get ganges path based on hostname
function [path2ganges] = get_ganges_path()

    path2ganges = '';

    [~, hostname] = system('hostname');

    if ~isempty(strfind(hostname, 'sallyw'))
        path2ganges = '/Volumes/';
    end

    if ~isempty(strfind(hostname, 'darya'))
        path2ganges = '/media/ganges/';
    end

    if ~isempty(strfind(hostname, 'spielwiese'))
        path2ganges = '/home/johannes/ganges/';
    end

    if strcmpi(path2ganges, '')
        here = pwd;
        i_ganges = strfind(here, '/ganges/');
        path2ganges    =  here(1:(i_ganges+6));
    end

end
