function interaction_score = calc_interaction_term(seg,broads,amp_or_del,thresh)
  
  %% Function calculates an interaction term which represents the
  %-logarithm of the odds ratio for having a focal segment of amplitude
  %af on top of a broad level of amplitude ab vs. on top of a broad level
  %of 0

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
  %amp_or_del > 0  for amplifications, < 0 for deletions
   
  if ~exist('thresh','var') || isempty(thresh)
    thresh = 0.1;
  end
    
  if amp_or_del > 0 
    interaction_score = -log(1);
  else
    broad_level = mean(full(broads(seg(2):seg(3),seg(5)))); %! approx. full roundoff w/SegArray
%!  broad_level = mean(broads(seg(2):seg(3),seg(5)));
    if broad_level >= thresh
      interaction_score = -log(max(0.1,1-broad_level-seg(12)));
    else
      interaction_score = 0;
    end
  end
