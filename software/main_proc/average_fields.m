function [avg] = average_fields(data, dt)
%%   [avg] = average_fields(data)
%        
%        This function averages all fields in 
%        data on time interval dt 
%
%
%        INPUT
%           avg.time  = matlab time 
%           dt        = time interval in [s] 
%
%
%   remodeled: 
%        Johannes Becherer
%        Wed Mar 28 15:41:54 PDT 2018


if ~isfield(data, 'time')
   if isfield(data, 'datenum')
      data.time = data.datenum(1:2:end);
   else
      error('there must be a time vector in the structure data');
      return;
   end
end

if ~isfield(data, 'time_tp')
   % for gust Tp is on the same time step as the other variables 
   if isfield(data, 'datenum')
      data.time_tp = data.datenum;
   else
      data.time_tp = data.time;
   end
end

if ~isfield(data, 'time_cmp')
   % usually there should always be a compass time, however 
   data.time_cmp = data.time;
end

%_____________________create new time vector______________________
sec_d = 3600*24; % factor sec to day

%_____________________time steps______________________
dt_time  =  median(diff(data.time));
dt_tp    =  median(diff(data.time_tp));
dt_cmp   =  median(diff(data.time_cmp));

%_____________________length______________________
N_time   =  length(data.time);
N_tp     =  length(data.time_tp);
N_cmp    =  length(data.time_cmp);

%_____________________average_window_width______________________
WW_time      = round(dt/(dt_time*sec_d));  
   if WW_time < 1; WW_time = 1; end;
   WW_time_half = round(WW_time/2);

WW_tp        = round(dt/(dt_tp*sec_d));  
   if WW_tp < 1; WW_tp = 1; end;
   WW_tp_half = round(WW_tp/2);

WW_cmp       = round(dt/(dt_cmp*sec_d));  
   if WW_cmp < 1; WW_cmp = 1; end;
   WW_cmp_half = round(WW_cmp/2);



%_____________________get all fields______________________
fs = fields(data);
   %---------------------initialize fileds fields----------------------
for f = 1:length(fs)
 N_field = length( data.(fs{f}) ); % field length
   if N_field ~= 1 % check if vector
      % ignor extra time vectors
      if ~(  strcmp(fs{f}, 'time_cmp') |  strcmp(fs{f}, 'time_tp'))

         if( strcmp(fs{f}, 'compass') | strcmp(fs{f}, 'cmp') | ...
             strcmp(fs{f}, 'pitch') | strcmp(fs{f}, 'roll') )
            % compass stuff is in deg and must be treated differently
            tmp         = data.(fs{f})/180*pi; % transfor into rad
            tmp         = exp(1i * tmp);  % convert to complex plain
            tmp         =  movmean( tmp, WW_cmp, 'omitnan');
            tmp         =  tmp( (1+WW_cmp_half):WW_cmp:(end-WW_cmp_half) );
            tmp         =  angle(tmp)*180/pi;
            tmp(tmp>360)=  tmp(tmp>360) - 360;
            tmp(tmp<0)  =  tmp(tmp<0)   + 360;
            avg.(fs{f}) = tmp;  

         else % normal fields
            if round(N_field/10) == round(N_time/10) % approx length time
               tmp         =  movmean( data.(fs{f}), WW_time, 'omitnan');
               avg.(fs{f}) =  tmp( (1+WW_time_half):WW_time:(end-WW_time_half) );
            elseif round(N_field/10) == round(N_tp/10) % approx length of time_tp
               tmp         =  movmean( data.(fs{f}), WW_tp, 'omitnan');
               avg.(fs{f}) =  tmp( (1+WW_tp_half):WW_tp:(end-WW_tp_half) );
            else % if there is no mtching time vector
               avg.(fs{f}) = data.(fs{f});
            end
         end
       

      end
   else % if scalar or structure
      avg.(fs{f}) = data.(fs{f});
   end
end  % for loop


end % function
