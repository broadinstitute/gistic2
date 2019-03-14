function [q,p,focal_d,ads,score_thresh,gg_rg,gg_ds] = compute_stats(D,scores,rg,score_type,alpha,qv_thresh,do_gene_gistic)
%COMPUTE_STATS convert marker scores into statistical P and Q values 
%
%  [Q,P,FOCAL_D,ADS,SCORE_THRESH,GG_RG,GG_DS] = compute_stats(D,...
%                     SCORES,RG,SCORE_TYPE,ALPHA,QV_THRESH,DO_GENE_GISTIC)
%
% INPUT ARGUMENTS:
%   D - copy number D_struct used by gene gistic to calculate gene scores
%   SCORES - a two-element amp/del cell array with each cell containing 
%     markers-by-samples arrays representing the (log) probability of each 
%     focal event under a background model of SCNA calculated from the 
%     aggregate of events.
%   RG - reference genome array of structs
%   SCORE_TYPE - a struct of scoring parameters: the 'res' field specifies
%     the CN resolution used in constructing empirical distributions
%   ALPHA - modeling parameters for exponential dependence of frequency on
%     amplitude
%   QV_THRESH - maximum threshold for significant q-values
%   DO_GENE_GISTIC - boolean set to 1 to do deletion analysis at gene
%     rather than marker level.
%
% OUTPUT ARGUMENTS:
%   Q - amp/del cell array of marker- or gene-level FDR "q values". The
%     dels cell will be a gene-level vector if the DO_GENE_GISTIC flag is
%     set.
%   P - amp/del cell array of marker- or gene-level P-values. The dels cell
%   will be a gene-level vector if the DO_GENE_GISTIC flag is set.
%   FOCAL_D - two-element cell array of the distribution of scores expected
%     by chance
%   ADS - two-element cell array of amplification (ads{1}) and deletion
%     (ads{2}) scores
%   SCORE_THRESH - two-element amp/del vector with the lowest score whose
%     q-value excedes QV_THRESH
%   GG_RG - if DO_GENE_GISTIC is set, this is a version of the rg reference
%     genome input to which additional fields have been added: 'snps' are
%     the marker indices for snps representing the gene, 'sc' is the gene
%     score, 'p' and 'q' are significance statistics fore the gene scores.
%     If do_gene_gistic is not set, then this output is empty
%   GG_DS - if DO_GENE_GISTIC is set, this is a cell array of the score
%     distributions for each gene size (in snp units). Each 
%     distribution is a normalized vector of probabilities for binned
%     scores. The width of each score bin is set by the score_type.res
%     input.
%     If DO_GENE_GISTIC is not set, then this output is empty

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


  verbose('Performing background permutations...',20);

  % Amps
  verbose('Amps...',30);
  [q{1},p{1},focal_d{1},ads{1},min_valA] = ...
      zigg_score_permutations(scores.amp+scores.aod,score_type.res*alpha(1));

  ads{1} = nanmean(scores.amp-min_valA,2);
  t = flipud(cumsum(flipud(focal_d{1})));
  p{1} = t(min(1+floor(ads{1}/(score_type.res*alpha(1))), ...
                     length(t)));
  q{1} = calc_fdr_value(p{1});
  
  % Dels
  verbose('Dels...',30);
  
  [q{2},p{2},focal_d{2},ads{2},min_valD] = ...
      zigg_score_permutations(scores.del+scores.doa,score_type.res*alpha(2));
  
  ads{2} = nanmean(scores.del-min_valD,2);
  t = flipud(cumsum(flipud(focal_d{2})));
  p{2} = t(min(1+floor(ads{2}/(score_type.res*alpha(2))), ...
                     length(t)));
  q{2} = calc_fdr_value(p{2});
  
  if do_gene_gistic
    verbose('Performing gene gistic...',20);
    
    [gg_rg, gg_ds] = gene_score_permutations(D, scores.del-min_valD,rg,alpha(2)* score_type.res);
    if isempty(gg_rg)
        do_gene_gistic = false;
    else
        gg_rg = order_rg_by_pos(gg_rg);

        gg_rg = gg_rg([gg_rg.p] >= 0);
        gene_p = [gg_rg.p];
        gene_q = calc_fdr_value(gene_p);
        gene_scores = score_from_p(gg_ds{1},gene_p)*score_type.res*alpha(2);
        for i=1:length(gg_rg)
          gg_rg(i).gene_scores = gene_scores(i);
          gg_rg(i).q = gene_q(i);
        end
        rr = find(gene_q<=qv_thresh);
        if ~isempty(rr)
          gene_score_thresh = min(gene_scores(rr));
        else
          gene_score_thresh = Inf;
          warning(['No significant Gene-GISTIC regions found!']); 
        end
    end
  else
    gg_rg = [];
    gg_ds = [];
  end

  % calculate marker-space score thresholds
  for k=1:2
    rr = find(q{k} <= qv_thresh);
    if ~isempty(rr)
      score_thresh(k)=min(ads{k}(rr));
    else
      score_thresh(k) = Inf;
      ampdel = {'amp','del'};
      warning(['No significant ' ampdel{k} ' regions found!']); 
    end
  end
  
  if do_gene_gistic
    score_thresh(2) = gene_score_thresh;
  end
  
  
  
