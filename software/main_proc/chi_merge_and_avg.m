function []  = chi_merge_and_avg(basedir, ddir, aw , time_lim, sname)
%% []  = chi_merge_and_avg(basedir, dir, aw, [time_lim, sdir])
%     
%     This function averages all idividual files in dir
%     in to a single file called 
%     dir_aw.mat
%
%     INPUT
%        basedir  :  base directory of unit
%        dir      :  sub directory of proc (e.g. temp, chi ...)
%        aw       :  average width in sec (if 0 no averaging)
%        time_lim :  optional, time limits to cut the field
%        sname     :  optional, path+file prefix for merged file
%
%   created by: 
%        Johannes Becherer
%        Tue Sep 20 16:28:51 PDT 2016
%

% if no time limits are set we use pratical no limits
if nargin < 4 | isempty(time_lim)
   time_lim = [datenum(1900,1,1) datenum(2100,1,1)];
end

if nargin < 5
    sname = [basedir 'proc' filesep ddir];
end

if nargin == 5
    sname = [basedir 'proc' filesep sname];
end

if ~isdir([basedir 'proc' filesep ddir])
   error([ddir ' is not a directory in ' basedir])
end

%---------------------find mat files----------------------
d = dir([basedir 'proc' filesep ddir]);
cnt = 1;
for i = 1:length(d)
   if ~(d(i).isdir)
      if d(i).name([-3:0]+end) == '.mat'
         fids{cnt} = d(i).name;
         cnt = cnt+1;
      end
   end
end

if cnt == 1  
   error([ 'There are no mat files in ' basedir ddir])
end


%---------------------first file----------------------
A = load([basedir 'proc' filesep ddir filesep fids{1}]);

   % what is the name of the field
   F  = fields(A);

%_____________________loop through substructures______________________
for s =1:length(F)

      % names of subfields
      FF = fields(A.(F{s}));

      avged = cell(size(fids));

      % average first file
      if aw == 0
         avged{1} = A.(F{s});
      else
         avged{1} = average_fields(A.(F{s}), aw);
      end

   %_____________________loop through rest of files______________________
   % for each file, read the data and append to cell array 'avged'
   for i = 2:length(fids)

      fid = fids{i};
      A = load([basedir 'proc' filesep ddir filesep fid]);

      % average
      if aw == 0
         avged{i} = A.(F{s});
      else
         avged{i} = average_fields(A.(F{s}), aw);
      end
   end

   % merge all fields
   avg.(F{s}) = merge_cell_structs(avged);
end

%---------------------save data----------------------
   if aw == 0
      sfid =[sname '.mat'] ;
   else
      sfid =[sname '_' num2str(aw) 'sec.mat'] ;
   end

   % find fields in avg
   Favg = fields(avg);
   for i =  1:length(Favg)
       if isfield(avg, 'time')
           eval([char(Favg(i)) ' =  time_lim_fields(avg.' char(Favg(i)) ', time_lim);']);
       else
           eval([char(Favg(i)) ' = avg.' char(Favg(i)) ';']);
       end
   end

   % save all fields in avg
   if length(Favg)==1
      save( sfid, char(Favg(1)), '-v7.3'); 
   elseif length(Favg)==2
      save( sffid, char(Favg(1)), char(Favg(2)), '-v7.3'); 
   elseif length(Favg)==3
      save( sfid, char(Favg(1)), char(Favg(2)), char(Favg(3)), '-v7.3'); 
   end

