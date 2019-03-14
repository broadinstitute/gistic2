function [nv,mv,stv]=dna_norm(v,flag)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

nc=size(v,2);
mv=nanmean(v,2);
stv=nanstd(v,0,2);
zero_idx=find(stv==0);		     
if ~isempty(zero_idx)
  stv(zero_idx)=1; % if std is zero the vector is const and nv
		   % will be zero
  warning('Zero std in dna_norm');
end
nv=(v-repmat(mv,1,nc))./(repmat(stv,1,nc));

if nargin>1 && flag
  nv=nv./sqrt(nc-1);
end
