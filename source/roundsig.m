function y=roundsig(x,dig)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if nargin==1
  dig=2;
end

if nnz(x~=0)==0
  y=zeros(size(x));
  return
end

s=(x<0);
x=abs(x);

scale=10.^(floor(log10(x))-dig+1);
y=scale.*round(x./scale);
y=y.*((-1).^s);
