function gene_ids = find_genes_from_snps(rg,snps)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)
  

  gene_ids = [];
  
  for i=1:length(rg)
    if ~isempty(intersect(rg(i).snps,snps))
      gene_ids = [gene_ids i];
    end
  end
  
  
