function y=interp_pwl(origx,f,x)
%INTERP_PWL Interpolate value of a function
%
%   Y = INTERP_PWL(ORIGX,F,X)
%
% Return interpolated Y value for X, using the vectors of (x,y) pairs in
% (ORIGX,F).
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

y=zeros(length(x),1);
for i=1:length(x)
  p_left=max(find(origx<x(i)));
  if isempty(p_left)
    p_left=1;
  end
  p_right=min(length(origx),p_left+1);
  if p_left == p_right || f(p_left) == f(p_right)
    y(i) = f(p_left);
  else
    r=(x(i)-origx(p_left)+eps)/(origx(p_right)-origx(p_left)+eps);
    y(i)=f(p_left)+r*(f(p_right)-f(p_left));
  end
end
