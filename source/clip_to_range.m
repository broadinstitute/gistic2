function x=clip_to_range(x,r,m);


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('m','var')
  x(x<r(1))=r(1);
  x(x>r(2))=r(2);
else
  x(x<r(1))=tanh((x(x<r(1))-r(1))/m(1))*m(1)+r(1);
  x(x>r(2))=tanh((x(x>r(2))-r(2))/m(2))*m(2)+r(2);
end

%if (x>r(2))
%  x=r(2);
%elseif x<r(1)
%  x=r(1);
%end
