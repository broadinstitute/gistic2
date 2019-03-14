function ads = smooth_gene_scores(D,gene_gistic_rg,gene_ads,ch)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  verbose('Smoothing score...',30);
  
  in_chr = find(D.chrn == ch);
  chr_zero = min(in_chr)-1;
  chr_max = max(in_chr);
  
  gene_gistic_rg = gene_gistic_rg(find([gene_gistic_rg.chrn] == ch));
  
  cur_sc = NaN(1,length(in_chr));
  
  for j=1:length(gene_gistic_rg)
    gene_snps = unique(gene_gistic_rg(j).snps - ...
                       chr_zero);
    cur_sc(gene_snps) = gene_ads(j);
  end
  
  no_gene = find(isnan(cur_sc));
  bpts = [1 find(diff(no_gene) > 1)]; 
  % if there is any intergenic space marked by NaNs
  if length(bpts) > 1
      % loop over intergenic spaces
      for j=1:length(bpts)
        % create index for begin/end/middle nogene spaces
        if j==1
          cur_snps = no_gene(bpts(j):bpts(j+1));
        elseif j==length(bpts)
          cur_snps = no_gene(bpts(j)+1);
        else
          cur_snps = no_gene(bpts(j)+1:bpts(j+1));
        end
        % assign maximum values for beginning/end/middle
        if min(cur_snps) == 1
          cur_sc(cur_snps) = cur_sc(max(cur_snps)+1);
        elseif max(cur_snps) == length(in_chr)
          cur_sc(cur_snps) = cur_sc(min(cur_snps)-1);
        else
          cur_sc(cur_snps) = max([cur_sc(min(cur_snps)-1) ...
                              cur_sc(max(cur_snps)+1)]);
        end
      end
  end
  
  ads = cur_sc;
  
