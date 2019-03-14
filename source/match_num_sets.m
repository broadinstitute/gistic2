function [Mt,m1,m2]=match_num_sets(set1,set2)
%MATCH_NUM_SETS index matching elements of two numerical vectors
%
% [MT,M1,M2] = match_num_sets(SET1,SET2)
%

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena.
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

Mt=[];
[uset1,~,uj1] = lunique(set1);
[uset2,~,uj2] = lunique(set2);

% check if both are unique
if length(uset1)==length(set1) && length(uset2)==length(set2) 
  
  c=intersect(set1,set2);
  
  m2=find(ismember(set2,c));
  x=find(ismember(set1,c));
  
  [~,s1i]=sort(set2(m2));
  [~,rs1i]=sort(s1i);
  
  [~,s2i]=sort(set1(x));
  [~,rs2i]=sort(s2i);
  
  m1=x(s2i(rs1i));
  Mt=sparse(m1,m2,true,length(set1),length(set2));
else
  [~,um1,um2]=match_num_sets(uset1,uset2);
  uMt=sparse(um1,um2,ones(length(um1),1),length(uset1),length(uset2));
  Mt=uMt(uj1,uj2);
  [m1,m2,~]=find(Mt);
end
