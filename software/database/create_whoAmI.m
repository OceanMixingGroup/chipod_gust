function [] = create_whoAmI( fname,  path2database, db_index)
%% [] = create_whoAmI( fname,  path2ganges, db_index)
%  
%     This function generate the whoAmI.m files
%
%     INPUT
%        fname       :  path to whoAmI.m  (e.g. '../mfiles/whoAmI.m');
%        path2ganges :
%        db_index    :  index of the instrument in the instruments database table
%
%
%   created by: 
%        Johannes Becherer
%        Sat Mar 24 17:16:18 PDT 2018



%_____________________get data from chosen entry______________________
% all fields in instrument table
I.fields{1}    = 'id';
I.fields{2}    = 'name';
I.fields{3}    = 'type';
I.fields{4}    = 'id_platform';
I.fields{5}    = 'start';
I.fields{6}    = 'stop';
I.fields{7}    = 'depth';
I.fields{8}    = 'datapath';
I.fields{9}    = 'comment';
I.fields{10}   = 'status';
I.fields{11}   = 'Owner_id';
I.fields{12}   = 'mab';

     % open db
     omg_db = sqlite(path2database);
        % instruments table
       for i = 1:length(I.fields)
           try 
               db_data = fetch(omg_db, ['select ' I.fields{i}  ' from instruments where id = ' num2str(db_index)] );
            catch me
               if contains(me.message, 'NULL')  % if null error
                  db_data = cell(1,1);
                  db_data{1,1} = '';
               else
                  error(me.message);
               end

            end
           if isnumeric((db_data{1,1}))
            I.values{i}   =  num2str(db_data{1,1});
           else
            I.values{i}   =  char(db_data{1,1});
           end
       end
     % close database
     close(omg_db);




% chipod or gust fields
if contains(I.values{3}, 'hipod')
   table = 'chipods';
   cnt = 1; 
   X.fields{cnt}  = 'sn'; cnt=cnt+1;
   X.fields{cnt}  = 'inst_id'; cnt=cnt+1;
   X.fields{cnt}  = 'depth'; cnt=cnt+1;
   X.fields{cnt}  = 'start'; cnt=cnt+1;
   X.fields{cnt}  = 'stop'; cnt=cnt+1;
   X.fields{cnt}  = 'problems'; cnt=cnt+1;
   X.fields{cnt}  = 'T1_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'T2_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'Tp1_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'Tp2_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'P_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'pitot_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'acc_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'cmp_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'dtime'; cnt=cnt+1;
   X.fields{cnt}  = 'dT1'; cnt=cnt+1;
   X.fields{cnt}  = 'dT2'; cnt=cnt+1;
   X.fields{cnt}  = 'dP'; cnt=cnt+1;
   X.fields{cnt}  = 'status'; cnt=cnt+1;
   X.fields{cnt}  = 'T1_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'T2_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'P_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'pitot_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'acc_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'cmp_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'T1_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'T2_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'T1P_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'T2P_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'P_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'Ax_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'Ay_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'Az_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'cmp_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'id'; cnt=cnt+1;

else     % gust
   table = 'gust';
   cnt = 1;
   X.fields{cnt}  = 'sn'; cnt=cnt+1;
   X.fields{cnt}  = 'type'; cnt=cnt+1;
   X.fields{cnt}  = 'start'; cnt=cnt+1;
   X.fields{cnt}  = 'stop'; cnt=cnt+1;
   X.fields{cnt}  = 'problems'; cnt=cnt+1;
   X.fields{cnt}  = 'T_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'Tp_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'P_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'pitot_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'acc_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'cmp_stop'; cnt=cnt+1;
   X.fields{cnt}  = 'dtime'; cnt=cnt+1;
   X.fields{cnt}  = 'dtemp'; cnt=cnt+1;
   X.fields{cnt}  = 'dP'; cnt=cnt+1;
   X.fields{cnt}  = 'inst_id'; cnt=cnt+1;
   X.fields{cnt}  = 'status'; cnt=cnt+1;
   X.fields{cnt}  = 'batteryfailure'; cnt=cnt+1;
   X.fields{cnt}  = 'todo'; cnt=cnt+1;
   X.fields{cnt}  = 'T_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'P_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'pitot_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'acc_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'cmp_sn'; cnt=cnt+1;
   X.fields{cnt}  = 'T_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'Tp_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'P_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'Ax_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'Ay_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'Az_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'cmp_cal'; cnt=cnt+1;
   X.fields{cnt}  = 'id'; cnt=cnt+1;
end

     % open db
     omg_db = sqlite(path2database);
        % instruments table
       for i = 1:length(X.fields)
           try 
               db_data = fetch(omg_db, ['select ' X.fields{i}  ' from ' table ' where inst_id = ' num2str(db_index)] );
            catch me
               if contains(me.message, 'NULL')  % if null error
                  db_data = cell(1,1);
                  db_data{1,1} = '';
               else
                  error(me.message);
               end

            end
           if isnumeric((db_data{1,1}))
            X.values{i}   =  num2str(db_data{1,1});
           else
            X.values{i}   =  char(db_data{1,1});
           end
       end
     % close database
     close(omg_db);

%_____________________write to whoAmI.m______________________

if exist(fname);
   overWrite =  input([ fname ' already exists. Shall I overwrite? y/n '] ,'s');
   if ~strcmp(overWrite, 'y')
      disp(' the file is not overwritten. action aborted!! ')
      return;
   end
end


fileID = fopen(fname,'w');
   fprintf(fileID, 'function [DB] = whoAmI() \n');
   fprintf(fileID, '%%  List of DB parameters for this particular unit\n');
   fprintf(fileID, '%%      automatically generated by whoAmI_generate.m\n\n\n');

   for i = 1:length(I.fields)
       line2add = ['  DB.instruments.' I.fields{i}];
       line2add(length(line2add)+1:32) = ' ';
       line2add   =  [line2add  '=  ''' I.values{i} ''';\n'];
       fprintf(fileID, line2add);
   end
   fprintf(fileID, '\n\n');
   for i = 1:length(X.fields)
       line2add = ['  DB.' table '.' X.fields{i}];
       line2add(length(line2add)+1:32) = ' ';
       line2add   =  [line2add  '=  ''' X.values{i} ''';\n'];
       fprintf(fileID, line2add);
   end
   fprintf(fileID, '\n\n');
   fprintf(fileID, 'end\n');
fclose(fileID);

end

