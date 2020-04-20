function [xd]=deglitch(x,npts,nstd,side, quiet_flag)
% function [xd]=deglitch(x,npts,nstd,side)
% x is the series to be deglitched
% npts - number of points to deglitch at a time 
%      - this shold be tuned the data
% nstd - number of std. dev. to remove outliers
% side - outliers will be removed only on positive 
%     (if side='+') or negative (if side='-')
%      side of the mean. If no side is defined, 
%      than outlier will be removed on both sides
% $Revision: 2.0 $  $Date: 2017/07/14 18:00:00 $
% Originally J. Moum

if nargin<3
    nstd=2;
elseif nargin<4
    side='b';
    quiet_flag = 0;
elseif nargin<5
    quiet_flag = 0;
end

np=npts; % no. of points to deglitch at a time
len=length(x);
xd=nan*ones(1,len);
nt=floor(len/np);

num_nan = sum(isnan(x));

try
    xd = x(1:np*nt);

    if size(xd, 2) ~= 1
        xd = xd';
    end

    sd = movstd(xd, npts, 'omitnan', 'endpoints', 'discard');
    mn = movmean(xd, npts, 'omitnan', 'endpoints', 'discard');

    % decimate to duplicate results of older version
    sd = sd(1:npts:end);
    mn = mn(1:npts:end);

    % make a 2d matrix and then flatten it to get mean, std in the
    % right places.
    sdmat = bsxfun(@times, sd, ones(1, npts))';
    sdvec = sdmat(1:np*nt)';
    mnmat = bsxfun(@times, mn, ones(1, npts))';
    mnvec = mnmat(1:np*nt)';

    if side == '+' | side == 'b'
        xd(xd > (mnvec + nstd*sdvec)) = NaN;
    end

    if side == '-' | side == 'b'
        xd(xd < (mnvec - nstd*sdvec)) = NaN;
    end

    % save memory
    clear mnmat sdmat mnvec sdvec

catch ME

    for it=1:nt
        xs=x((np*(it-1)+1):np*it);
        idnan=isnan(xs);
        xm=nanmean(xs);
        xsd=nanstd(xs);
        if side=='b'
            %        id=find((xs)>(xm+(nstd*xsd)) | (xs)<(xm-(nstd*xsd)));
            id=find(abs(xs-xm)>nstd*xsd);
        elseif side=='+'
            id=find(xs>(xm+(nstd*xsd)));
        elseif side=='-'
            id=find(xs<(xm-(nstd*xsd)));
        else
            disp('Wrong "side" is defined')
            break
        end
        xs(id)=NaN;
        xd((np*(it-1)+1):np*it)=xs;
    end
end

% now do the tail end points
xs=x(np*nt+1:len);
xm=nanmean(xs);
xsd=nanstd(xs);
if side=='b'
    id=find(abs(xs-xm)>nstd*xsd);
elseif side=='+'
    id=find(xs>(xm+(nstd*xsd)));
elseif side=='-'
    id=find(xs>(xm-(nstd*xsd)));
else
    disp('Wrong "side" is defined')
end
xs(id)=NaN;
xd(np*nt+1:len)=xs;
if size(xd,1)~=size(x,1); xd=xd'; end

if ~quiet_flag
    num_newnan = sum(isnan(xd));
    disp(['deglitch removed ' num2str(num_newnan-num_nan) ' points = ' ...
          num2str((num_newnan-num_nan)/(length(x) - num_nan) * 100, '%.1f') ' % | ' ...
          'avg(NaNed points) = ' num2str(nanmean(x(isnan(xd))), '%.1e') ' | ' ...
          'median(NaNed points) = ' num2str(nanmedian(x(isnan(xd))), '%.1e')])
end