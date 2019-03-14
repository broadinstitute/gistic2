function x=find_text_pos(a,w,s,L,U)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

a=as_column(a);
w=w/2;
s=s/2;
n=length(a);

H=zeros(2*n,2*n);
for i=1:n
    H((2*i-1):(2*i),(2*i-1):(2*i))=0.5;
end
f=-[a a]';
f=f(:);

A=zeros(n-1,2*n);
for i=1:(n-1)
    A(i,2*i)=1;
    A(i,2*i+1)=-1;
end
b=-2*s*ones(n-1,1);

Aeq=zeros(n,2*n);
for i=1:n
    Aeq(i,2*i-1)=-1;
    Aeq(i,2*i)=1;
end
beq=2*w*ones(n,1);

LB=L-w*ones(2*n,1);
UB=U+w*ones(2*n,1);

X0=a(1)+(0:(n-1))'*(a(end)-a(1))/(n-1);
optimoptions = optimset('Display','off');  %suppress "optimization terminated" verbosity
x=quadprog(H,f,A,b,Aeq,beq,LB,UB,[],optimoptions);
x=reshape(x,2,length(x)/2);
x=mean(x)';
