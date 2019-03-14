function D=jaccard(M,N)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if nargin==1
  n=size(M,1);
  m1=double(M==1);
  en=m1*m1';
  sr=sum(m1,2);
  dn=repmat(sr,1,n)+repmat(sr',n,1)-en;
  
  m2=double(isnan(M));
  if nnz(m2)>0
    dn=dn-m1*m2'-m2*m1';
  end
  
  D=en./dn;
else
  m=size(M,1);
  n=size(N,1);
  m1=double(M==1);
  n1=double(N==1);
  
  en=m1*n1';
  srm=sum(m1,2);
  srn=sum(n1,2);
  dn=repmat(srm,1,n)+repmat(srn',m,1)-en;
  
  m2=double(isnan(M));
  n2=double(isnan(N));
  if (nnz(m2)>0) | (nnz(n2)>0)
    dn=dn-m1*n2'-m2*n1';
  end
  
  D=en./dn;
  
end



