function c=combine_cyto(cyto,level)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

s={cyto.name};
switch level
 case 'subcyto'
 case 'cyto'
  s=regexprep(s,'\..*','');
 case 'arm'
  s=regexprep(s,'[\.0-9]+$','');  
 case 'chr'
  s=regexprep(s,'[pq][\.0-9]+$','');  
end
[u,ui,uj]=unique_keepord(strvcat(s),'rows');
for i=1:length(ui)
  j=find(uj==i);
  c(i)=cyto(ui(i));
  c(i).start=cyto(min(j)).start;
  c(i).end=cyto(min(j)).end;
  c(i).name=deblank(u(i,:));
  c(i).stain='Empty';
end
