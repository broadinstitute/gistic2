function [regs range_dist] = calculate_wide_peaks_for_regs(Z,Qs,regs,perm_ads,conf_level,score_thresh,score_type,range_struct,gene_gistic_rg,peak_types)
    

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if ~exist('gene_gistic_rg','var') || isempty(gene_gistic_rg)
    do_gene_gistic = 0;
  else
    do_gene_gistic = 1;
  end
    
  if ~exist('peak_types','var') || isempty(peak_types)
    peak_types = {'robust'};
  end
  
  if ~isempty(strmatch('robust',peak_types))
    if ~exist('range_struct','var') || ~isstruct(range_struct)
      error(['must include structure specifying range distribution ' ...
             'parameters.']);
    else
      switch range_struct.method
       case 'permute'
        do_permutations = 1;
        range_dist = range_stuct.range_dist;
       case 'approximate'
        do_permutations = 0;
        max_chr_length = max(arrayfun(@(x) length(find(Z.chrn == x)),unique([Z.chrn])));
        max_win_length = min(max_chr_length,floor(length(Z.chrn)/2));
        win_range = 1:max_win_length;
        if isfield(range_struct,'interpl_method')
          interpl_method = range_struct.interpl_method;
        else
          interpl_method = 'linear';
        end
        approx_range = interpolate_range_dist_percentile(perm_ads,score_type,conf_level,win_range,interpl_method);
      end
    end
  end
  
  for k=1:2
    if k==1
      Q = Qs.amp;
    else
      Q = Qs.del;
    end
    if ~isempty(strmatch('robust',peak_types)) && do_permutations
      if ~exist('range_dist','var') || isempty(range_dist{k})
        range_dist{k} = cell(1,length(find(Z.chrn==1)));  
      end
    end
    
    if k==2 && do_gene_gistic
      gene_gistic = 1;
    else
      gene_gistic = 0;
    end
    
    for i=1:size(regs{k},2)
      temp_reg = regs{k}(i);
      cur_ads = regs{k}(i).ads;
      ch = regs{k}(i).chrn;
      in_chr = find(Z.chrn == ch);
      chr_zero = min(in_chr)-1;
      chr_max = max(in_chr);
      cur_score = temp_reg.score;
      temp_reg.peak_st = temp_reg.peak_st - chr_zero;
      temp_reg.peak_en = temp_reg.peak_en - chr_zero;

      if gene_gistic
        
        genes_in_chr = find([gene_gistic_rg.chrn] == ch);
        gene_st = temp_reg.robust_reg_st_gene;
        if gene_st == 0
          temp_reg.robust_reg_st = chr_zero+1;
          temp_reg.robust_reg_st_gene = 1;
        else
          temp_reg.robust_reg_st = ...
              min(gene_gistic_rg(genes_in_chr(gene_st)).snps);
        end
        gene_en = temp_reg.robust_reg_en_gene;
        if isinf(gene_en)
          temp_reg.robust_reg_en = chr_max;
          temp_reg.robust_reg_en_gene = length(cur_ads);
        else
          temp_reg.robust_reg_en = ...
              max(gene_gistic_rg(genes_in_chr(gene_en)).snps);
        end
           
        cur_ads = smooth_gene_scores(Z,gene_gistic_rg(genes_in_chr), ...
                                     cur_ads,ch);
      end
          
      if ~isempty(strmatch('robust',peak_types)) 
        if do_permutations
          [robust_wide_st robust_wide_en thresh range_dist{k}] = robust_peak(temp_reg,cur_ads,perm_ads{k},cur_score-score_thresh(k),score_type,temp_reg.robust_reg_st-chr_zero,temp_reg.robust_reg_en-chr_zero,conf_level,range_dist{k});
        else
          [robust_wide_st robust_wide_en thresh] = robust_peak(temp_reg,cur_ads,perm_ads{k},cur_score-score_thresh(k),score_type,temp_reg.robust_reg_st-chr_zero,temp_reg.robust_reg_en-chr_zero,conf_level,[],[],approx_range{k});
        end
        
        regs{k}(i).robust_peak_wide_st = max(chr_zero+1,robust_wide_st+chr_zero);
        regs{k}(i).robust_peak_wide_en = min(chr_max,robust_wide_en+chr_zero);
      end
      
      %% leave 1-out
      if ~isempty(strmatch('loo',peak_types))
        x = Z.dat{k}(in_chr,:);
        y = x>0;
        rl = runlength(y);
        samples = find(regs{k}(i).samples);
        rli = find_rl(rl(samples),[regs{k}(i).peak_st-min(in_chr)+1 regs{k}(i).peak_en-min(in_chr)+1]);
        
        z=zeros(size(x));
        for l=1:length(samples)
          for j=1:length(rli{l})
            rng=rl{samples(l)}(rli{l}(j),1):rl{samples(l)}(rli{l}(j),2);
            z(rng,samples(l))=x(rng,samples(l));
          end
        end
        [regs{k}(i).loo_wide_st,regs{k}(i).loo_wide_en]=add_wide_peak(z,chr_zero);
      end 
      
      if strmatch('robust',peak_types{1})
        regs{k}(i).peak_wide_st = regs{k}(i).robust_peak_wide_st;
        regs{k}(i).peak_wide_en = regs{k}(i).robust_peak_wide_en;
      else
        regs{k}(i).peak_wide_st = regs{k}(i).loo_wide_st;
        regs{k}(i).peak_wide_en = regs{k}(i).loo_wide_en;
      end
      
      if k==1
        disp(['Amplification wide peak at ' genomic_location(Z,{regs{k}(i).peak_wide_st:regs{k}(i).peak_wide_en})]);
      else
        disp(['Deletion wide peak at ' genomic_location(Z, ...
                                                        {regs{k}(i).peak_wide_st:regs{k}(i).peak_wide_en})]);
      end
      
    end                     
    
  end %% end of k=1:2 loop
