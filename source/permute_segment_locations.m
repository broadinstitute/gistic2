function perm_ads = permute_segment_locations(Qs,nperm,nsamples,nsnps,clean_struct)
% permute_segment_locations -- Permutes segment locations to
% create null scores for computation of range distribution
% 
%
% perm_ads = permute_segment_locations(Qs,nperm,nsamples,nsnps,clean_struct)
%    
%
%
% Qs (requird): 1x2 cell array containing ziggurat segments for amps/dels
% nperm(required) = number of permutations to perform for amps/dels
% nsamples(required) = number of samples in dataset
% nsnps(required) = number of snps in dataset
% clean_struct (optional): a structure which can be used to clean drivers
%                from data before permutations, if desired.  
%                If used, this structure must contain the following fields:
% 
%                          clean_up = 1 (use clean-up) or 0 (no clean-up)
%                          seg_assignment = 1x2 cell array containing assignment of each
%                                segment to significant peaks in data
%                          score_type = score_type struct from gistic run
%                          d = background distribution
%                          alpha = alpha used for scoring segments
%                          (default to [1 1])
% 
%  Outputs: perm_ads, a 1x2 cell array.  Each field is itself a 1 x nperm cell array 
%  of permuted background scores for amps/dels, respectively.
%
% ---
% $Id$
% $Date: 2014-01-31 15:30:25 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  
  
  if ~exist('clean_struct','var') || isempty(clean_struct)
    clean_up = 0;
  else
    clean_up = clean_struct(1).clean_up;
    if clean_up
      
      if isfield(clean_struct,'seg_assignment')
        seg_assignment{1} = clean_struct(1).seg_assignment;
        seg_assignment{2} = clean_struct(2).seg_assignment;
      else
        error(['Must have seg_assignment field in clean_struct to use clean_up']);
      end
      
      if isfield(clean_struct,'d')
        d{1} = clean_struct(1).d;
        d{2} = clean_struct(2).d;
      else
        error(['Must have background distribution in clean_struct to use ' ...
               'clean-up']);
      end
      
      if isfield(clean_struct,'score_type')
        score_type = clean_struct(1).score_type;
      else
        error('Must have score_type in clean_struct to use clean-up!');
      end
      
      if isfield(clean_struct,'alpha')
        alpha = clean_struct.alpha;
      else
        alpha = [1 1];
      end
    
    end
  end
    
  perm_ads = cell(1,2);
  
  for k=1:2
    if k==1
      Q = Qs.amp;
    else
      Q = Qs.del;
    end
    
    if clean_up
      %% Compute mean of bg_distribution 
      gmean = 0;
      for i=1:length(d{k})
        gmean = gmean + i*d{k}(i);
      end
      gmean = gmean*score_type.res*alpha(k);
      
    
      %% Find segments assigned to no significant peaks -- these
      %are all passengers
      
      nsegs = size(seg_assignment{k},1);
      un_assigned = find(sum(seg_assignment{k},2) == 0);
      assigned = setdiff(1:nsegs,un_assigned);
      
      %% For remaining assigned peaks, determine probability that each
      %is a passenger
    
      reg_disjoint_scores = sum(seg_assignment{k},1)/nsamples;
      weight_matrix = zeros(size(seg_assignment{k}));
      weight_matrix(assigned,:) = seg_assignment{k}(assigned,:) > 0;
      
      prob_pass = ones(nsegs,1);
      prob_pass(assigned) = gmean./(weight_matrix(assigned,:)*reg_disjoint_scores');
    else
      prob_pass = ones(nsegs,1);
    end
      
    new_ads = cell(1,nperm);
    
    for i=1:nperm
      disp(['Performing permutation ' num2str(i) ' of ' num2str(nperm)]);
      coin = rand(nsegs,1);
      keep_segs = find(coin <= prob_pass);
      cur_ads = permute_seg_locations(Q(keep_segs,:),nsamples,nsnps);
      new_ads{i} = cur_ads;
    end

    perm_ads{k} = new_ads;
  
  end
  
