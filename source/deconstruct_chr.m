function [ZA ZD] = deconstruct_chr(B,chr_bpt,p_level,q_level)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if ~exist('p_level','var') || isempty(p_level)
    p_level = 0;
  end
  
  if ~exist('q_level','var') || isempty(q_level)
    q_level = 0;
  end
  
  if ~exist('chr_bpt','var') || isempty(chr_bpt)
    chr_bpt = max(B(:,3));
  end
  
  bpt_row = find(B(:,3) == chr_bpt);
  BP = [];
  BQ = [];
  if isempty(bpt_row)
    error('Chromosome breakpoint must correspond to segment breakpoint!');
  elseif bpt_row == size(B,1)
    BP = B(1:bpt_row,:);
    BP(:,4) = BP(:,4) - p_level;
  elseif bpt_row == 1
    BQ = B(bpt_row:end,:);
    BQ(:,4) = BQ(:,4) - q_level;
  else
    BP = B(1:bpt_row,:);
    BQ = B(bpt_row+1:end,:);
    BP(:,4) = BP(:,4) - p_level;
    BQ(:,4) = BQ(:,4) - q_level;
  end
  
  % prepare for zigg_deconstruction
  [BPa BPd] = prepare_B(BP);
  [BQa BQd] = prepare_B(BQ);
    
  % zigg_deconstruction on both p and q arms
  ZAp = add_broad_levels_to_zigg(atomic_zigg_deconstruction(BPa),p_level);
  ZDp = add_broad_levels_to_zigg(atomic_zigg_deconstruction(BPd),p_level);
  ZAq = add_broad_levels_to_zigg(atomic_zigg_deconstruction(BQa),q_level);
  ZDq = add_broad_levels_to_zigg(atomic_zigg_deconstruction(BQd),q_level);
  
  % merge ziggurats
  ZA = cat(1,ZAp,ZAq);
  ZD = cat(1,ZDp,ZDq);
