function [ZA,ZD] = deconstruct_sample(B,broad_levels,chr_bpts)
  
%% Loop over chromosomes

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  ZAt = cell(1,max(B(:,1)));
  ZDt = ZAt;
  
  for ch=1:max(B(:,1))
    cur_segs = find(B(:,1) == ch);
    Bt = B(cur_segs,:);
    [ZAt{ch} ZDt{ch}] = deconstruct_chr(Bt,chr_bpts(ch),broad_levels(2* ...
                                                      ch-1),broad_levels(2*ch));
  end
      
  ZA = cat(1,ZAt{:});
  ZD = cat(1,ZDt{:});
  
