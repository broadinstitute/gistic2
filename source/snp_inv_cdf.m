function g=snp_inv_cdf(d,score_type,conf_level)
% Function takes g-score distribution d, score_type variable, and the desired percentile (conf_level), and
% returns the g-score associated with the desired percentile (g{1} for amp, g{2} for del)
% ---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  
  if isempty(d)
    error('Must supply snp_score background distribution!')
  end
  
  if isempty(score_type)
    error('Must supply score_type variable!')
  end
  
  if isempty(conf_level) || ~isfloat(conf_level) || conf_level < 0 || ...
        conf_level > 1
    error('Must supply a confidence level between 0 and 1!')
  end
  
  %% Convert from probability distribution d to cumulative distribution t
  if size(d,2)>1
    for k = 1:size(d,2)
      t{k} = flipud(cumsum(flipud(d{k}))); 
  
  %% Now find the index associated with probability level 
  
      ind{k} = min(find(t{k} <= 1-conf_level))-1;
      
      if isempty(ind{k})
        ind = length(t{k});
      end
      
    
  %% Finally, find the g-score associated with that index level
    
      g{k} = ind{k}*score_type.res;
    
    end
  else 
    t = flipud(cumsum(flipud(d)));
    ind = min(find(t<=1-conf_level))-1;
    if isempty(ind)
      ind=length(t);
    end
    g = ind*score_type.res;
  end
  
