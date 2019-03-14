function hd=hist2d(x,y,rx,ry)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

xn=zeros(size(x));
yn=zeros(size(x));
xn(find(x<=rx(1)))=1;
yn(find(y<=ry(1)))=1;
xn(find(x>=rx(end)))=length(rx);
yn(find(y>=ry(end)))=length(ry);

for i=1:(length(rx)-1)
  xn(find( x >=rx(i) & x <rx(i+1)))=i;
end
for i=1:(length(ry)-1)
  yn(find( y >=ry(i) & y <ry(i+1)))=i;
end

%hd=crosstab_full(xn,yn,1:max(length(rx)+1,length(ry)+1));
hd=crosstab_full(xn,yn,1:max(length(rx),length(ry)));
hd=hd(1:length(rx),1:length(ry));
