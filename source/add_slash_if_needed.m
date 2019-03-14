function dirst=add_slash_if_needed(dirst)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if isempty(dirst)
  dirst=filesep;
else
  if dirst(end)~=filesep
    dirst=[dirst filesep];
  end
end
