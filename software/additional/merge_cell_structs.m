% [merged] = merge_cell_structs(in, concat_2d)
%           in : cell array of structures
%    concat_2d : if 0 (default), this routine combines 1d vectors to make a long
%                1d vector
%                if 1, this routine combines 1d vectors to yield a 2d array
function [merged] = merge_cell_structs(in, concat_2d)

    if ~exist('concat_2d', 'var'), concat_2d = 0; end

    mat = cell2mat(in);
    F = fields(mat);

    for f=1:length(F)
        if isstruct(in{1}.(F{f}))
            merged.(F{f}) = merge_struct_array([mat.(F{f})]);
        else
            [~, concat_dim] = max(size(mat(1).(F{f})));
            if concat_2d
                if concat_dim == 1
                    concat_dim = 2;
                else
                    concat_dim = 1;
                end
            end
            merged.(F{f}) = cat(concat_dim, mat.(F{f}));
        end
    end
end