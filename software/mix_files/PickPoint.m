function [date] = PickPoint(n)

    if ~exist('n', 'var'), n = 1; end
    [x, y] = ginput(n);
    datestr(x)
