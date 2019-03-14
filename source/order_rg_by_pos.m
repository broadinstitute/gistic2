function rg=order_rg_by_pos(rg)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~isfield(rg,'chrn')
  rg = add_chrn(rg);
end

tmp=cat(1,rg.chrn)*1e11+0.5*(double(cat(1,rg.start))+double(cat(1,rg.end)));
[stmp,si]=sort(tmp);

rg=rg(si);
