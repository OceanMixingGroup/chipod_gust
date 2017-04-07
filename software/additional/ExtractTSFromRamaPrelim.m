function [T1, T2] = ExtractTSFromRamaPrelim(ramaname, ChipodDepth)

    rama = load(ramaname);
    keyboard;

    [z1,z2] = ChooseDepthLevels(rama.depth, ChipodDepth)

    T1.time = rama.time;
    T1.S = rama.sal(z1,:);
    T1.T = rama.temp(z1,:);
    T1.z = rama.depth(z1);

    T2.time = rama.time;
    T2.S = rama.sal(z2,:);
    T2.T = rama.temp(z2,:);
    T2.z = rama.depth(z2);

end

% find instruments with least amount of missing data within 15m of chipod
% missing data checked over chipod deployment time
function [index1, index2] = ChooseDepthLevels(depth, ChipodDepth)

    index1 = find( (abs(-depth + ChipodDepth) < 7.5) ...
                    & (ChipodDepth > depth));

    index2 = find( (abs(-depth + ChipodDepth) < 7.5) ...
                    & (ChipodDepth < depth));

    assert(depth(index1) < ChipodDepth & depth(index2) > ChipodDepth, ...
           'Measurements don''t span Chipod depth!');
end
