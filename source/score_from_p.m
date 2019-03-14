function scores=score_from_p(d,p)
% Takes a distribution d and vector of p-values p to calculate scores
% that would have generated each p-value
%
%---
% $Id$
% $Date: 2014-01-31 15:30:44 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


scores = zeros(length(p),1);
  
t = flipud(cumsum(flipud(d)));

for i = 1:length(p)

  scores(i) = min([find(t<p(i));length(t)+1]);

end

