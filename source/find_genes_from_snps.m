function gene_ids = find_genes_from_snps(rg,snps)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  gene_ids = [];
  
  for i=1:length(rg)
    if ~isempty(intersect(rg(i).snps,snps))
      gene_ids = [gene_ids i];
    end
  end
  
  
