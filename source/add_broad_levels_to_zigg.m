function Z = add_broad_levels_to_zigg(Z,broad_level)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if ~isempty(Z)
    Z(:,6) = Z(:,6) + broad_level;
    Z(:,7) = Z(:,7) + broad_level;
  end
