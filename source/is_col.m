function res=is_col(rc)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

switch rc(1:min(length(rc),3))
 case {'col','sam','con','exp'}
  res=1;
 case {'row','gen','mir','mar','snp'}
  res=0;
 otherwise
  error('no match');
end
