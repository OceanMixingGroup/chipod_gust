function [merged] = merge_cell_structs(in)

    mat = cell2mat(in);
    F = fields(mat);

    for f=1:length(F)
        if isstruct(in{1}.(F{f}))
            merged.(F{f}) = merge_struct_array([mat.(F{f})]);
        else
            [~,i_max] = max(size(mat(1).(F{f})));
            merged.(F{f}) = cat( i_max, mat.(F{f}));
        end
    end
end