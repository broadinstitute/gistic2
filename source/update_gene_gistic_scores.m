function gene_scores = update_gene_gistic_scores(Z,gene_gistic_rg, ...
                                                 gene_gistic_ds,res,chr_zero)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if ~exist('chr_zero','var') || isempty(chr_zero)
    chr_zero = 0;
  end
  
  verbose('Updating gene gistic scores for %d genes',20, ...
          length(gene_gistic_rg));
  
  gene_scores = zeros(1,length(gene_gistic_rg));
  
  for i=1:length(gene_gistic_rg)
    temp_score = max(0,mean(max(Z(gene_gistic_rg(i).snps-chr_zero,:),[],1)));
    nsnps = gene_gistic_rg(i).nsnps;
    temp_t = flipud(cumsum(flipud(gene_gistic_ds{nsnps})));
    temp_p = temp_t(min(1+floor(temp_score/res),length(temp_t)));
    gene_scores(i) = score_from_p(gene_gistic_ds{1},temp_p)*res;
  end
