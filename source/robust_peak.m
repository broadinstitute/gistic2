function [wide_peak_start wide_peak_end thresh range_dist gene_list] = robust_peak(reg,ads,perm_ads,score_thresh,score_type,region_start,region_end,conf_level,range_dist,gene_gistic_rg,approx_range)
% robust_peak -- Compute wide peak for a region using robust
% 
%
% [temp_start temp_end thresh range_dist] = robust_peak(reg,ads,perm_ads,conf_level,score_thresh,score_type,region_start,region_end,range_dist)
%
%
% reg(required): current reg structure for which wide_peak will be added
%                At minimum, this must contain peak start, peak stop, and
%                peak score
% ads (required): The current score profile on which robust should be run
% perm_ads (required): A cell array containing random permutations of the
%                      gistic score (used for computing
%                      range_distribution)
% score_thresh (required): The threshold score associated with the
%                          region_start/region_end, it marks the region
%                          that would have independently peeled off
%                          as significant in the data
% score_type (required): score_type structure used in GISTIC
% region_start, region_end (required): the boundaries of the region in
%                                      which to search for thewide_peak, usually
%                                      defined by the limits of score_thresh
%   
% conf_level (default = 0.75): The level of confidence for assigning wide
%                             peak
% range_dist (required): A cell array containing the range distribution for each
%             window size
%
% Outputs: wide_peak_start/wide_peak_end = boundaries of robust wide peak
%          thresh = score at boundary (where difference between peak and
%                                      it is consistent with range_distribution)
%          range_dist = returns the range_distribution with (potentially)
%                       new window size computed (speed-up so same range dist doesn't
%                       have to be computed each time)
% ---
% $Id$
% $Date: 2014-01-31 15:30:37 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  
  if exist('approx_range','var') & ~isempty(approx_range)
    do_permutations = 0;
  else
    do_permutations = 1;
  end
  
  if ~exist('gene_gistic_rg','var') || isempty(gene_gistic_rg)
    gene_gistic = 0;
  else
    gene_gistic = 1;
  end
  
  if ~exist('conf_level','var') || isempty(conf_level)
    conf_level = 0.75;
  end
  
  peak_score = reg.score;
  
  peak_start = reg.peak_st;
  peak_end = reg.peak_en;
  num_peak_snps = Inf;
  temp_start = region_start;
  temp_end = region_end;
  num_temp_snps =  temp_end-temp_start+1;
  num_perms = 0;
  
  if gene_gistic
    peak_gene_start = min(reg.peak_genes);
    region_gene_start = reg.robust_reg_st_gene;
    peak_gene_end = max(reg.peak_genes);
    region_gene_end = reg.robust_reg_en_gene;
    temp_gene_start = region_gene_start;
    temp_gene_end = region_gene_end;
  end
  
  while num_temp_snps < num_peak_snps
    num_perms = num_perms+1;
    num_peak_snps = num_temp_snps;
       
    if do_permutations
      if isempty(range_dist{num_temp_snps})
        verbose('Permutation #%d with window size of %d',20,[num_perms num_temp_snps]);
        range_dist{num_temp_snps} = calc_range_dist_in_windows(perm_ads,num_temp_snps,score_type.res);
      else
        verbose('Permutation #%d: using old permuted range distribution for window size of %d',20,[num_perms,num_temp_snps]);
      end
      gt_index = snp_inv_cdf(range_dist{num_temp_snps}',score_type,conf_level);
    else
      gt_index = approx_range(num_temp_snps);
    end
      
    thresh = peak_score-gt_index;
    
    if thresh <= score_thresh  %% Wide peak is all of 'region'
      temp_start = region_start;
      temp_end = region_end;
      num_temp_snps = region_end-region_start+1;
      num_peak_snps = num_temp_snps;
      thresh = score_thresh;
    else %% Find wide peak region and update number of temporary snps
      if ~gene_gistic
        rg = intersect(region_start:region_end,find(ads < thresh));
        temp_start = max(intersect(rg,region_start:peak_start));
        if isempty(temp_start)
          temp_start = region_start;
        end
        temp_end = min(intersect(rg,peak_end:region_end));
        if isempty(temp_end)
          temp_end = region_end;
        end
      else
        rg_gene = intersect(region_gene_start:region_gene_end,find(ads ...
                                                          < thresh));
        temp_gene_start = max(intersect(rg_gene,region_gene_start: ...
                                        peak_gene_start));
        if isempty(temp_gene_start)
          temp_gene_start = region_gene_start;
        end
        temp_gene_end = min(intersect(rg_gene,peak_gene_end: ...
                                      region_gene_end));
        if isempty(temp_gene_end)
          temp_gene_end = region_gene_end;
        end
        temp_start = min(gene_gistic_rg(temp_gene_start).snps);
        temp_end = max(gene_gistic_rg(temp_gene_end).snps);
      end
      num_temp_snps = temp_end-temp_start+1;
    end
  end
    
  wide_peak_start = temp_start;
  wide_peak_end = temp_end;

  if gene_gistic
    gene_list = unique({gene_gistic_rg(temp_gene_start: ...
                                       temp_gene_end).symb});
  end
