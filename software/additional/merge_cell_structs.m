function [merged] = merge_cell_structs(in)

    mat = cell2mat(in);
    F = fields(mat);

    for f=1:length(F)
        if isstruct(in{1}.(F{f}))
            merged.(F{f}) = merge_struct_array([mat.(F{f})]);
        else
            merged.(F{f}) = [mat.(F{f})];
        end
    end
end