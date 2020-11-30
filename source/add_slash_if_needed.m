function dirst=add_slash_if_needed(dirst)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)


if isempty(dirst)
  dirst=filesep;
else
  if dirst(end)~=filesep
    dirst=[dirst filesep];
  end
end
