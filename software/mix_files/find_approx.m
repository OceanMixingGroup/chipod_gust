function g = find_approx(m,v,n)
%     g = find_approx(m,v,[n])
% find index (g) of matrix (m) that is most nearly equal to value (v); this
% is similar to g=find(m==v), except that the nearest approximate equality
% is found if no exact equality exists.
% the third argument (default n=1) tells how many values to find; n=3 means
% the nearest 3 indices in order of descending nearness.

if nargin>2
    g=nan(n,1);
    for nn=1:n
        [nul g(nn)]=min(abs(m-v));
        m(g(nn))=nan;
    end

else
    [nul g]=min(abs(m-v));
    if isnan(nul), g=nan; end
end