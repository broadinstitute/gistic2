function str=fill_empty(str,estr)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if iscell(str)
  e=find(cellfun('isempty',str));
  str(e)={estr};
else
  if isempty(str)
    str=estr;
  end
end
