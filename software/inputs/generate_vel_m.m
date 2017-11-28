function [vel_m] = generate_vel_m(time, z_adcp, u, v, z_chi,  sdir)
%% [vel_m] = generate_vel_m(time, z_adcp, u, v, z_chi, sdir)
% 
%        This function generates an input file for chi-processing vel_m.mat
%        directory sdir
%        
%        Input
%           time        : adcp time
%           z_adcp      : adcp depth 
%           u           : east velocity matrix
%           v           : north velocity matrix
%           z_chi       : instrument depth 
%           sdir        : directory to save vel_m.mat to (if empty no saving)
% 
%        Output
%           vel_m.time  : time vector  
%           vel_m.u     : east velocity vector
%           vel_m.v     : north velocity vector
%           vel_m.z_chi : instrument depth 
%
%   created by: 
%        Johannes Becherer
%        Fri Sep  2 11:20:31 PDT 2016


% init variables
vel_m.u     = nan(size(time));
vel_m.v     = nan(size(time));
vel_m.spd     = nan(size(time));
vel_m.time  = time;
vel_m.depth = z_chi;


if (length(z_chi) == 1) & (z_adcp == z_chi)   % in case it is already interpolated

   vel_m.u = u;
   vel_m.v = v;
   
else     % in cse it is at a differnt depth interpolate
   % determin the size of z_adcp (matrix? or vector?)
   sz = size(z_adcp);

   % if vector generate matrix
   if sz(1)~=length(time) & sz(2)~=length(time)
      [z_adcp, ~] = meshgrid(z_adcp, time);
   end

   if sz(2) == length(time) % fix dimension order  
       z_adcp = z_adcp';
   end
   if size(u,2) == length(time)
       u = u';
       v = v';
   end

   if length(z_chi) == 1
       z_chi = z_chi * ones(size(time));
   end

   % loop through every time step
   for t = 1:length(time)
      z_tmp = z_adcp(t,:);
      [z_sort, ii_sort_z] = sort(z_tmp); 

      % look if z_chi is out of bound
      if z_chi(t) < z_sort(1)
         vel_m.u(t) = u(t, ii_sort_z(1));
         vel_m.v(t) = v(t, ii_sort_z(1));
      elseif z_chi(t) > z_sort(end)
         vel_m.u(t) = u(t, ii_sort_z(end));
         vel_m.v(t) = v(t, ii_sort_z(end));
      else % if not out of z bound
         vel_m.u(t) = interp1(z_tmp, u(t,:), z_chi(t));
         vel_m.v(t) = interp1(z_tmp, v(t,:), z_chi(t));
      end

   end
end

vel_m.spd = hypot(vel_m.u, vel_m.v);
vel_m.U = vel_m.u + 1i * vel_m.v;
vel_m.comment = '(u,v) = velocities; spd = speed; U = u+iv';

if ~isempty(sdir)
   save([sdir 'vel_m.mat'], 'vel_m');
end
