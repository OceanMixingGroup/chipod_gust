function [IDstr] = fix_underscore( ID )
% [IDstr] = fix_underscore( ID )
%
%       ID = some sort of string that had underscores
%       IDstr = the same string as ID except with '_' replaced with '\_'
%
%   When printing strings with underscores, matlab creates subscripts. 
%   This function allows underscores to print as underscores.


nn = 0;
for tt = 1:length(ID)
    nn = nn+1;
    if strcmp(ID(tt),'_')
        IDstr(nn) = '\';
        IDstr(nn+1) = '_';
        nn = nn+1;
    else
        IDstr(nn) = ID(tt);
    end
end


end

