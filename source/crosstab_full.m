function ct=crosstab_full(c1,c2,v)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

n=length(v);
ct=zeros(n,n);
for i=1:length(v)
  if isnan(v(i))
    fc1(isnan(c1))=i;
    fc2(isnan(c2))=i;    
  else
    fc1(c1==v(i))=i;
    fc2(c2==v(i))=i;
  end
end

for i=1:length(fc1)
  ct(fc1(i),fc2(i))=ct(fc1(i),fc2(i))+1;
end



