function s=add_chrn(s)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~isfield(s,'chrn')
  if length(s)==1
    s.chrn=chromosome2num(s.chr);
  else
    for i=1:length(s)
      s(i).chrn=chromosome2num(s(i).chr);
    end
  end
end
