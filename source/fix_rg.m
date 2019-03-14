function [rg,mi] = fix_rg(rg,mi)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if size(rg,1) > size(rg,2)
    rg = rg';
  end
  
  if any(diff(rg) > 1)
     warning(['the peak has more than one segment']);
     segs = [0 find(diff(rg)~=1) length(rg)];
     [sx si] = max(diff([0 find(diff(rg)>1) length(rg)]));
     rg = rg(segs(si)+1:segs(si+1));  
     mi=rg(1);
  end
  
  
