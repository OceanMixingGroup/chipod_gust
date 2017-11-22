function [out]=calibrate_polynomial(in,coef)
% function [out]=calibrate_poly(in)
%   in - input data structure
%   calibrates using polynomial coefficients in head.coef
%   $Revision: 1.1.1.1 $  $Date: 2011/05/17 18:08:05 $

% in case pavan forgot to put the whole 5 coef vector
if numel(coef) == 2
    coef = [coef 0 0 0];
end
    

out=coef(1)+coef(2).*in+coef(3).*(in.^2)+coef(4).*(in.^3);