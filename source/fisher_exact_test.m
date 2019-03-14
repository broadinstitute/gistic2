function [p1,p2]=fisher_exact_test(a,b,c,d,eps,tail_st)
% [p1,p2]=fisher_exact_test(a,b,c,d,eps,tail_st)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~exist('eps','var') || (exist('eps','var') && isempty(eps))
    eps=1e-6;
end

if exist('tail_st','var') && strcmp(lower(tail_st),'left')
  left_side=1;
else
  left_side=0;
end

if ~exist('b','var') || (exist('b','var') && isempty(b))
  b=a(1,2);
  c=a(2,1);
  d=a(2,2);
  a=a(1,1);
end
  
if size(a,1)>1
  for i=1:size(a,1)
    [p1(i),p2(i)]=fisher_exact_test(a(i),b(i),c(i),d(i),eps);
  end
else    
  p1=0;  
  if ~left_side
    if (a+b)==0 || (a/(a+b)>c/(c+d))
      ta=a;
      tb=b;
      a=c;
      b=d;
      c=ta;
      d=tb;
    end
  end

  const=gammaln(a+b+1)+gammaln(c+d+1)+gammaln(a+c+1)+gammaln(b+d+1)-gammaln(a+b+c+d+1);
  for i=a:-1:max(0,a-d)    
    %    delta=exp(chooseln(i,a+b)+chooseln(c+a-i,c+d)-chooseln(a+b,a+b+c+d));
    delta=exp(const-gammaln(i+1)-gammaln(b+a-i+1)-gammaln(c+a-i+1)-gammaln(d-a+i+1));
    %    delta
    p1=p1+delta;
    if (delta/p1)<eps
      break
    end
  end
  %p1

  if ~left_side
    p2=p1;
    
    j=ceil(((a+b)*(a+2*c)-a*(c+d))/(a+b+c+d));
    
    for i=j:1:min(a+c,a+b)
      %    delta=exp(chooseln(i,a+b)+chooseln(c+a-i,c+d)-chooseln(a+b,a+b+c+d));
      delta=exp(const-gammaln(i+1)-gammaln(b+a-i+1)-gammaln(c+a-i+1)- ...
                gammaln(d-a+i+1));
      %    delta
      p2=p2+delta;
      if (delta/p2)<eps
        break
      end
    end
    % p2
  else
    p2=NaN;
  end
end


