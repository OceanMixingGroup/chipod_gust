% [merged] = merge_struct_array(in, concat_2d)
%           in : array of structures
%    concat_2d : if 0 (default), this routine combines 1d vectors to make a long
%                1d vector
%                if 1, this routine combines 1d vectors to yield a 2d array
function [merged] = merge_struct_array(in, concat_2d)

    if ~exist('concat_2d', 'var'), concat_2d = 0; end

    FF = fields(in);
    for ff=1:length(FF)
        [~, concat_dim] = max(size(in(1).(FF{ff})));
        if concat_2d
            if concat_dim == 1
                concat_dim = 2;
            else
                concat_dim = 1;
            end
        end

        merged.(FF{ff}) = cat(concat_dim, in.(FF{ff}));
        %merged.(FF{ff}) = [in.(FF{ff})];
    end
end