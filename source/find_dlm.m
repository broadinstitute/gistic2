function d=find_dlm(fname,dlm)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('dlm','var') || isempty(dlm)
  dlm=[ char(9) ',|'];
end

f=fopen(fname,'r');
l=fgetl(f);
for i=1:length(dlm)
  h(i)=length(find(l==dlm(i)));
end
[hm,hi]=max(h);
d=dlm(hi);
