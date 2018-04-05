clear all;
close all;


% find ganges on computer (This assumes you operate on ganges)
% if you are not on ganges set path2ganges manually
path2ganges = get_ganges_path;
path2database  =  [path2ganges '/work/database/omg.sqlite'];

addpath(genpath('./chipod_gust/software/'));

%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

   %unit = '52';

%_____________________Ask the data base______________________
% open omg database
omg_db = sqlite(path2database);

% find all instruments in database that belong to mooring
sql_string = ['SELECT inst_name, inst_id,  datapath, inst_type, platform_name ' ...
                     'FROM instOnPlat where inst_name like "%' unit '%" and ' ...
                     '( inst_type like "%chipod%" or '...
                     'inst_type like "%gust%"  )'];
   %disp('The following SQL statement is applied to the omg_database:');
   %disp(['  ' sql_string]);
data = fetch(omg_db, sql_string );

% close database
close(omg_db);

%_____________________loop through all matches______________________
N_inst   = size(data,1) ;
for i = 1:N_inst
   inst_name{i}      = char(data{i,1});
   inst_id(i)        = data{i,2};
   datapath{i}       = char(data{i,3});
   inst_type{i}      = char(data{i,4});
   platform_name{i}  = char(data{i,5});
end


%_____________________choose an instrument______________________
if ~isempty(N_inst)
   disp('   I found the following matches in the database');
   for i =1:N_inst
      disp(['      ' num2str(i) '     : sn = ' inst_name{i} ' ; platform = ' platform_name{i} ' ; path = ' datapath{i} ]);
   end
else
   disp('   I did not find any matching instrument ');
   disp('    options:'); 
end
disp('      0     : create a new entry in the database');
disp('      9999  : cancel')

prompt   = ['Which option do you choose : ']; 
inkey    =  input(prompt);


%_____________________what to do______________________
% new entry
if inkey == 0  
   disp(' A new entry in the database is created');
   disp('_____________________');
   
   newSN    =  input('What is the instruments Serial number? ' ,'s');
   InstTyps{1} =  'chipod';
   InstTyps{2} =  'gust (straight)';
   InstTyps{3} =  'gust (gooseneck)';
   InstTyps{4} =  'gust (AoA)';
   InstTyps{5} =  'gust (MultiGusT)';
   disp('What type of instrument?');
   for i = 1:length(InstTyps)
      disp(['   ' num2str(i) '  : ' InstTyps{i}]);
   end
   
   type_id = 0;
   while ~sum(type_id == [1:5]) % ask until a valid choise is done
      type_id  = input('choose between 1-5 : ');
   end

   i_ganges =  strfind(basedir, '/ganges/');
   datapath =  basedir(30:end);


   % open db
   omg_db = sqlite(path2database);


      % create new entry in instruments table
      genStr_inst = ['insert into instruments (name, id_platform,  type, datapath) values ("' newSN '", 138, "' InstTyps{type_id} '","' datapath  '")'];
	   exec(omg_db, genStr_inst);

      %get the last entry in instruments table
      db_data = fetch(omg_db, 'select max(id) from instruments' );
      db_index   =  db_data{1,1};

      % write correponding entry in chipods/gust table
      if type_id == 1 % chipods
         genStr_chi = ['insert into chipods (sn, inst_id) values ("' newSN '","' num2str(db_index) '")'];
	      exec(omg_db, genStr_chi);
      else     % gust
         genStr_gust = ['insert into gust (sn, type,  inst_id) values ("' newSN  '","' InstTyps{type_id} '","' num2str(db_index) '")'];
	      exec(omg_db, genStr_gust);
      end

   % close database
   close(omg_db);

   disp('A new instrument was written to the database!!')

% choose existing entry
elseif inkey <= N_inst  
   db_index = inst_id(inkey);
   disp(['You chose option No: ' num2str(inkey)])
      i = inkey;
      disp(['      ' num2str(i) '     : sn = ' inst_name{i} ' ; platform = ' platform_name{i} ' ; path = ' datapath{i} ]);


% cancel action
else   
   disp('You cancled the generation of the whoAmI.csv file');
   return;
end


create_whoAmI( [basedir '/mfiles/whoAmI.m'],  path2database, db_index);

disp(' whoAmI.m sucessfully generated !')

