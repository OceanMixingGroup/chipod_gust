function [hash] = githash(fname, gitdir)

    if ~exist('fname', 'var')
        fname = '.';
    else
        fname = ['-- ' fname];
    end

    if ~exist('gitdir', 'var')
        gitdir = '--git-dir=./chipod_gust/.git/';
    else
        gitdir = ['--git-dir=' gitdir];
        if isempty(findstr(gitdir, '.git'))
            gitdir = [gitdir '/.git/'];
        end
    end

    [badexit, hashout] = system(['TERM=xterm git ' gitdir ...
                        ' log -n 1 --no-color --pretty=format:''%H'' ' fname ' < /dev/null']);
    if badexit
        % needed on MATLAB servers :(
        [~, hashout] = system(['git ' gitdir ...
                        ' log -n 1 --no-color --pretty=format:''%H'' ' fname ' < /dev/null']);
        hash = strtrim(hashout(7:46));
        return
    end

    try
        % remove bash escape characters
        hash = strtrim(hashout(8:47));
    catch ME
        warning('githash failed');
    end
