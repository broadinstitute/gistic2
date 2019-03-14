function rli=find_rl(rl,pos,above_val)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('above_val','var') || isempty(above_val)
  above_val=-Inf;
end
if iscell(rl)
  for i=1:length(rl)
    rli{i}=find_rl(rl{i},pos);
  end
else
  if length(pos)==1
    rli=find(pos>=rl(:,1) & pos<=rl(:,2) & rl(:,3)>above_val);
  else
%    rli=find(((rl(:,1)>=pos(1) & rl(:,1)<=pos(2)) | (rl(:,2)>=pos(1) & rl(:,2)<=pos(2))) & (rl(:,3)>above_val));
    rli=find(((rl(:,2)>=pos(1) & rl(:,1)<=pos(2))) & (rl(:,3)>above_val));
  end   
end
