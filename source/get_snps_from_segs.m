function snps = get_snps_from_segs(Q)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  snps = [];
  
  for i=1:size(Q,1)
    snps = unique([snps Q(i,2):Q(i,3)]);
  end
  
  
