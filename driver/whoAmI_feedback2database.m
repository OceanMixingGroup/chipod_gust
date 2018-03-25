clear all;
close all;

% find ganges on computer (This assumes you operate on ganges)
%   if you are not on ganges set path2ganges manually
here = pwd;
i_ganges = strfind(here, '/ganges/');
path2ganges =  here(1:i_ganges);


addpath(genpath('./chipod_gust/software/'));

%_____________________get info from local file______________________
DB       = whoAmI();
inst_id  = DB.instruments.id;

%_____________________find difference______________________
create_whoAmI([pwd '/db_test.m'], path2ganges, inst_id );

%disp the following changes would be written to the data base
!diff whoAmI.m db_test.m
   overWrite =  input(['Are you sure that you want to submit these changes to the database? y/n '] ,'s');
   if ~strcmp(overWrite, 'y')
      disp(' Nothing was written to database. action aborted!! ')
      !rm db_test.m
      return;
   end
!rm db_test.m




%_____________________feed back______________________
tables = fields(DB);

% open db
omg_db = sqlite([path2ganges 'ganges/work/database/omg.sqlite']);

for t = 1:length(tables)
   columns  =  fields(DB.(tables{t}));
   for c = 1:length(columns)

      if ~strcmp(columns{c}, 'id') % don't update id 
			if strcmp(tables{t}, 'instruments') 
				sqlstr   =  ['UPDATE ' tables{t} ...
							 ' SET ' columns{c} ' = "' DB.(tables{t}).(columns{c})  '"' ...
							' WHERE id = ' inst_id ];
			else
				sqlstr   =  ['UPDATE ' tables{t} ...
							 ' SET ' columns{c} ' = "' DB.(tables{t}).(columns{c})  '"' ...
							' WHERE inst_id = ' inst_id ];
			end
	      exec(omg_db, sqlstr);
	   end

   end
end

% close database
close(omg_db);
