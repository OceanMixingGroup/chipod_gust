function [fids, fdate] = chi_find_rawfiles(basedir, t_lims)
%% [fids, fdate] = chi_find_rawfiles(basedir, t_lims)
%     This function finds all the raw data that match unit
%     in rawdir, and returns a list of files  (fids{i})
%     and the corresponding date stamps (as string, fdate{i})
   
if nargin < 2
   t_lims = [datenum(2000, 1, 1, 0, 0, 0) datenum(2060, 1, 1, 0, 0, 0)];
end


   rawdir = [basedir filesep 'raw' filesep];
   unit   = chi_get_unit_name(basedir);

   d = dir(rawdir);
   cnt = 1;
   fids = [];
   for i = 1:length(d)
      if(~d(i).isdir)
         if(d(i).name([1:3])=='raw') % check if raw-file
            if ~isempty(strfind(d(i).name, unit)) % check if correct unit
               % read in the file date from file name
               is1 = strfind(d(i).name,'_');  % find under score
               is2 = strfind(d(i).name,'.');  % find dot

               fdate{cnt} = d(i).name((is1+1):(is2-1));
               fids{cnt} = d(i).name;
               cnt = cnt+1;
            else
               disp([d(i).name ' does not fit unit ' unit]);
            end
         else % for some older deployments, no 'raw_' at beginning of filenames
             if ~strcmp(d(i).name(1),'.')
                 if(d(i).name([-2:0]+end)== unit([-2:0]+end)) % check if correct unit
                    % read in the file date from file name
                    is1 = 0;  % no underscore
                    is2 = strfind(d(i).name,'.');  % find dot

                    fdate{cnt} = d(i).name((is1+1):(is2-1));
                    fids{cnt} = d(i).name;
                    cnt = cnt+1;
                 else
                    disp([d(i).name ' does not fit unit ' unit]);
                 end  
             end
         end
      end
   end



       %  only days within time limits
         dateformat = 'yymmddHHMM';
         for d = 1:length(fdate)
            ftime(d) = datenum( fdate{d}, dateformat(1:length(fdate{d})));
         end
         dd = find( ftime>=floor(t_lims(1)) & ftime <= ceil(t_lims(2)) );
         if length(dd) < length(fids)
            for d = 1:length(dd)
               tmpf{d} = fids{dd(d)};
               tmpd{d} = fdate{dd(d)};
            end
            fids  = tmpf;
            fdate = tmpd;
         end


   if(isempty(fids))
      error(['There are no matching raw files at the given path ' rawdir])
      return;
   end
   
end   
