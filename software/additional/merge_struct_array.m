function [merged] = merge_struct_array(in)

    FF = fields(in);
    for ff=1:length(FF)
        [~,i_max] = max(size(in(1).(FF{ff})));
        merged.(FF{ff}) = cat( i_max, mat.(FF{ff}));
        %merged.(FF{ff}) = [in.(FF{ff})];
    end
end