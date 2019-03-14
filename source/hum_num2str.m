function st=hum_num2str(num,use_commas)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

st=num2str(num);
if use_commas
  str=[];
  for i=(length(st)-2):-3:2
    str=[st(i:(i+2)) ',' str];
  end
  str=[st(1:(i-1)) ',' str];
  str=str(1:(end-1));
  st=str;
end

