function [data] = apply_time_lims_raw_data(data, TL)

   if isfield(data, 'T1') % chhpod
       data.T1 = apply_limits(data.time, data.T1, TL.T1);
       data.T2 = apply_limits(data.time, data.T2, TL.T2);

       data.T1Pt = apply_limits(data.time_tp, data.T1Pt, TL.Tp1);
       data.T2Pt = apply_limits(data.time_tp, data.T2Pt, TL.Tp2);
   else % gust
       data.T = apply_limits(data.time, data.T, TL.T);

       data.TPt = apply_limits(data.time_tp, data.TPt, TL.Tp);
   end

    data.W = apply_limits(data.time, data.W, TL.pitot);

    data.P = apply_limits(data.time, data.P, TL.P);
    data.depth = apply_limits(data.time, data.depth, TL.P);
    data.p_dis_z = apply_limits(data.time, data.p_dis_z, TL.P);
    data.p_vel_z = apply_limits(data.time, data.p_vel_z, TL.P);

    data.AX = apply_limits(data.time, data.AX, TL.acc);
    data.AY = apply_limits(data.time, data.AY, TL.acc);
    data.AZ = apply_limits(data.time, data.AZ, TL.acc);
    data.AXtilt = apply_limits(data.time, data.AXtilt, TL.acc);
    data.AYtilt = apply_limits(data.time, data.AYtilt, TL.acc);
    data.AZtilt = apply_limits(data.time, data.AZtilt, TL.acc);
    data.Acc = apply_limits(data.time, data.Acc, TL.acc);
    data.a_dis_x = apply_limits(data.time, data.a_dis_x, TL.acc);
    data.a_vel_x = apply_limits(data.time, data.a_vel_x, TL.acc);
    data.a_dis_y = apply_limits(data.time, data.a_dis_y, TL.acc);
    data.a_vel_y = apply_limits(data.time, data.a_vel_y, TL.acc);
    data.a_dis_z = apply_limits(data.time, data.a_dis_z, TL.acc);
    data.a_vel_z = apply_limits(data.time, data.a_vel_z, TL.acc);

    data.cmp = apply_limits(data.time_cmp, data.cmp, TL.cmp);
end

function [in] = apply_limits(in_time, in, time_lim)

    ind1 = find_approx(in_time, time_lim(1));
    ind2 = find_approx(in_time, time_lim(2));

    if ~isempty(ind1), in(1:ind1) = NaN; end
    if ~isempty(ind2), in(ind2:end) = NaN; end
end
