function [scores D params] = score_genome(D,focals,params)
%SCORE_GENOME - transform SCNA event lengths and amplitudes into scores
%
%    [SCORES D PARAMS] = score_genome(D,FOCALS,PARAMS)
%
% The returned SCORES is a two-element amp/del cell array with each cell
% containing markers-by-samples arrays representing the (log)
% probability of each focal event under a background model of SCNA
% calculated from the aggregate of events. Currently, each segment
% amplitude is multiplied by a scaling factor. For deletions, an
% interaction term is added that represents the conditional probability of
% observing a focal deletion on top of a broad deletion. 
%
% The input D is a copy number structure which must contain the Qs field of
% ziggurat segments as well as copy number data. The returned D contains a
% modified Qs field with the SCNA scores in column 12 calculated from the
% amplitudes in column 12 of each event set.
%
% FOCALS is a structure containing copy number data matrices for
% amplification (.amp), deletion (.del), amplification over
% deletion (.aod) and deletion over amplification (.doa).
%
% PARAMS is a structure containing optional parameter fields
%    PARAMS.t_amp             Amplification filtering threshold.
%    PARAMS.t_del             Deletion filtering threshold. 
%    PARAMS.broad_len_cutoff  Filtering parameter for broad events (arm-scaled)
%    PARAMS.alpha             Exponential scaling factors ([amp del]) for modeling
%                             event frequency as a function of length. If empty,
%                             alpha will be derived from the data.
%    PARAMS.pseudocount       Percentage of total amplification or deletion events 
%                             to be added to each event bin in the alpha calculation
%                             (if alpha is calculated)
% The returned PARAMS has had all unspecified fields filled in with their
% default values.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


  params = impose_default_value(params,'t_amp',0.1);
  params = impose_default_value(params,'t_del',0.1);
  params = impose_default_value(params,'broad_len_cutoff',0.98);
  params = impose_default_value(params,'alpha',[2.5145 2.1653]);
  params = impose_default_value(params,'use_segarray',false);
  params = impose_default_value(params,'pseudocount',0);
  
  Qs = D.Qs;
    
  if isempty(params.alpha)
    % Compute alphas from data
    params.alpha = [2.5145 2.1653];
  elseif length(params.alpha) == 1
    params.alpha = [params.alpha params.alpha];
  end
  
  % Score segments
  verbose('Scoring segments...',20);
  Qs.amp(:,4) = params.alpha(1)*Qs.amp(:,12);
  Qs.aod(:,4) = params.alpha(1)*Qs.aod(:,12);
  Qs.del(:,4) = params.alpha(2)*Qs.del(:,12);
  Qs.doa(:,4) = params.alpha(2)*Qs.doa(:,12);

  % Add broad interaction term (to deletions only)
  broads = reconstruct_genomes(rmfield(Qs,{'aod','doa'}),...
                               struct('broad_or_focal','broad',...
                                      'broad_len_cutoff',params.broad_len_cutoff,...
                                      't_amp',0,...
                                      't_del',0,...
                                      'column_to_add',12,...
                                      'use_segarray',params.use_segarray,...
                                      'rows',length(D.pos),...
                                      'cols',length(D.sdesc)));

  
  focal_del_segs = find(Qs.del(:,8) < params.broad_len_cutoff & Qs.del(:,12) >= params.t_del);
  focal_doa_segs = find(Qs.doa(:,8) < params.broad_len_cutoff & Qs.doa(:,12) >= params.t_del);
  
  %% Focals Dels
  for i=1:length(focal_del_segs)
    interaction_term = calc_interaction_term(Qs.del(focal_del_segs(i),:),broads.del,-1,params.t_del);
    Qs.del(focal_del_segs(i),4) = Qs.del(focal_del_segs(i),4)+interaction_term;
  end
  
  for i=1:length(focal_doa_segs)
    interaction_term = calc_interaction_term(Qs.doa(focal_doa_segs(i),:),broads.del,-1,params.t_del);
    Qs.doa(focal_doa_segs(i),4) = Qs.doa(focal_doa_segs(i),4) + interaction_term;
  end
  
  % Scoring genomes
  verbose('Scoring genomes...',20);
  scores = reconstruct_genomes(Qs,struct('broad_or_focal','focal',...
                                         'broad_len_cutoff',params.broad_len_cutoff,...
                                         't_amp',params.t_amp,...
                                         't_del',params.t_del,...
                                         'column_to_add',4,...
                                         'use_segarray',params.use_segarray,...
                                         'rows',length(D.pos),...
                                         'cols',length(D.sdesc)));
  
  % Subtract log(alpha) from scoresA,scoresD
  scores.amp = scores.amp-log(params.alpha(1));
  scores.del = scores.del-log(params.alpha(2));
  % return Qs in the D struct
  D.Qs = Qs;
  
