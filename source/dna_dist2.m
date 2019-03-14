function d=dna_dist2(vecs)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

N=size(vecs,1);
g2=sum(vecs.*vecs,2);
[U,S,V]=svd(vecs,0);
d=U*(U'.*repmat((-2)*diag(S).^2,1,size(U,1)));
st=1;
while (st < N)
  en=min(st+1000,N);
  d(st:en,:)=d(st:en,:)+repmat(g2(st:en),1,N)+repmat(g2',en-st+1,1);
  st=en+1;
end
d(d<0)=0;
st=1;
while (st < N)
  en=min(st+1000,N);
  d(st:en,:)=sqrt(d(st:en,:));
  st=en+1;
end

% d=repmat(g2,1,N)+d+repmat(g2',N,1);
