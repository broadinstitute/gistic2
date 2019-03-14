function [broad_level max_Q max_score num_levels] = find_max_broad_level_by_table(B,hd,xamp,ylen,arm_fract)
      

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  unique_levels = unique(B(:,4)); %% Add zero as a level?
  num_levels = length(unique_levels);
  
  if length(unique_levels)>1
    max_level = max(unique_levels);
    min_level = min(unique_levels);
    [QAS QDS] = ziggurat_on_extremes(B,min_level,max_level,hd,xamp,ylen);
    
    for k=1:length(unique_levels)
      Q{k} = iterative_ziggurat(QAS,QDS,unique_levels(k),hd,xamp,ylen);
      Qbroad = [B(1,1) min(B(:,2)) max(B(:,3)) unique_levels(k) B(1,5) 0 ...
                unique_levels(k) arm_fract 0];
      Qbroad(:,9) = score_ziggs_by_table(Qbroad,hd,xamp,ylen);
      
      score(k) = sum(Q{k}(:,9))+Qbroad(:,9);

      if Qbroad(:,4) ~= 0
        Q{k} = cat(1,Q{k},Qbroad);
      end
    end      
            
    [max_score mk] = max(score);
    broad_level = unique_levels(mk);  
    max_Q = Q{mk}; 
  elseif length(unique_levels) == 1
    broad_level = unique_levels(1);
    max_Q = [B(1,1:3) broad_level B(1,5) 0 broad_level arm_fract];
    max_score = score_ziggs_by_table(max_Q,hd,xamp,ylen);
  else
    broad_level = 0;
    max_Q = [];
    max_score = 0;
  end
