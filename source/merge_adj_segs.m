function M = merge_adj_segs(M)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)
  
  % Takes a dataarray M corresponding to chromosomal segments and merges
  % adjacent rows whose value in amp field are equal
  % If ids variable is provided, only checks rows immediately after rows
  % in ids.  Defaults to all rows if ids is not provided.

    i=1;
    while i < size(M,1)
      if M(i,4) == M(i+1,4) % check whether amps are equal
        M(i,3) = M(i+1,3); % change ending snp
        M(i,6) = M(i,6)+M(i+1,6); % add fractions
        M=M(setdiff(1:size(M,1),i+1),:); % remove redundant segment
      else
        i=i+1;
      end
    end
    
        
