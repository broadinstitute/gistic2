function snps = get_snps_from_segs(Q)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)
  

  snps = [];
  
  for i=1:size(Q,1)
    snps = unique([snps Q(i,2):Q(i,3)]);
  end
  
  
