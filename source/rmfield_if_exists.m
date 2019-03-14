function x=rmfield_if_exists(s,flds)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ischar(flds)
  if isfield(s,flds)
    x=rmfield(s,flds);
  else
    x=s;
  end
else
  x=s;
  for i=1:length(flds)
    x=rmfield_if_exists(x,flds{i});
  end
end
