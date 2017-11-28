function [S] = time_lim_fields(S, t_lims)
%%    [S] = time_lim_fields(S, t_lims)
%     
%     This function cuts all time dependend fielde in structure
%     S such that they fit in time limits set by t_lims
%
% 
%     Note this function assumes there is a field called time in the
%     the structure and that time dependend fields have max 4 demensions
%
%     created by
%        Johannes Becherer


% get all fields in S
fs = fields(S);

% what is the size of the time vector assuming their is a time vector
N   = length(S.time); 
% find all indexes that lie in the time limits
ii  = find( S.time >= t_lims(1) & S.time <= t_lims(2) );

% loop through all fields in S
for i = 1:length(fs)
   
   % find the demension that belongs to time

   d_time = find( size(S.(fs{i})) == N );
   if ~isempty(d_time)  % treat field only iof it has a time denesion
      switch d_time
         case 1
            switch ndims(S.(fs{i}))
               case 1
                  S.(fs{i}) = S.(fs{i})(ii);
               case 2
                  S.(fs{i}) = S.(fs{i})(ii,:);
               case 3
                  S.(fs{i}) = S.(fs{i})(ii,:,:);
               case 4
                  S.(fs{i}) = S.(fs{i})(ii,:,:,:);
             end
         case 2
            switch ndims(S.(fs{i}))
               case 2
                  S.(fs{i}) = S.(fs{i})(:,ii);
               case 3
                  S.(fs{i}) = S.(fs{i})(:,ii,:);
               case 4
                  S.(fs{i}) = S.(fs{i})(:,ii,:,:);
             end
         case 3
            switch ndims(S.(fs{i}))
               case 3
                  S.(fs{i}) = S.(fs{i})(:,:,ii);
               case 4
                  S.(fs{i}) = S.(fs{i})(:,:,ii,:);
             end
         case 4
            S.(fs{i}) = S.(fs{i})(:,:,:,ii);
      end

   end

end
