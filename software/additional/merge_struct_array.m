function [merged] = merge_struct_array(in)

    FF = fields(in);
    for ff=1:length(FF)
        merged.(FF{ff}) = [in.(FF{ff})];
    end
end