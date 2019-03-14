function ha = clean_hist(ha,num_snps)
%CLEAN_HIST remove lowest scoring snps from a histogram
%
%   HA = CLEAN_HIST(HA,NUM_SNPS)
%
% NUM_SNPS lowest scoring SNPS are removed from the histogram represented
% by the vector HA.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  bin = 1; % Current bin to remove from
  total = 0; % Total number of snps removed
  
  while total < num_snps
    if total+ha(bin) >= num_snps %% if there are enough snps in the
                                      %current bin
      ha(bin) = ha(bin)-(num_snps-total);
      total = num_snps;
    else %% there are not enough snps in the current bin!
      total = total+ha(bin);
      ha(bin) = 0;
      bin = bin+1;
    end
  end
      
      
