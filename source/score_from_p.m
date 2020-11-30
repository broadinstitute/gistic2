function scores=score_from_p(d,p)
% Takes a distribution d and vector of p-values p to calculate scores
% that would have generated each p-value

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

scores = zeros(length(p),1);
t = flipud(cumsum(flipud(d)));
for i = 1:length(p)
  scores(i) = min([find(t<p(i));length(t)+1]);
end

