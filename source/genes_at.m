function [idx,closest_idx]=genes_at(rg,chr,pos_start,pos_end,closest_flag,partial_hits)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

closest_idx=[];

if pos_end<pos_start
  error('End position is before start position');
end

if ~isfield(rg,'chrn')
  rg=add_chrn(rg);
end

in_chr=find(cat(1,rg(:).chrn)==chr);

if (exist('partial_hits','var') && (partial_hits==1 || isempty(partial_hits)) ...
    || ~exist('partial_hits','var')) 
  in_reg=find((cat(1,rg(in_chr).start)<=pos_end) & (cat(1,rg(in_chr).end)>= ...
                                                    pos_start));
else
  in_reg=find((cat(1,rg(in_chr).start)>=pos_start) & (cat(1,rg(in_chr).end)<= ...
                                                    pos_end));
end

idx=in_chr(in_reg);

if isempty(idx) && exist('closest_flag','var') && ~isempty(closest_flag)
  before_idx=find((cat(1,rg(in_chr).end)<pos_start));
  if ~isempty(before_idx)
    [before_dist,before_min_idx]=min(int32(pos_start)-cat(1,rg(in_chr(before_idx)).end));
  else
    before_dist=Inf;
  end

  after_idx=find((cat(1,rg(in_chr).start)>pos_end));
  if ~isempty(after_idx)
    [after_dist,after_min_idx]=min(cat(1,rg(in_chr(after_idx)).start)-int32(pos_end));
  else
    after_dist=Inf;
  end
  
  if before_dist < after_dist
    closest_idx=in_chr(before_idx(before_min_idx));
    closest_dist=before_dist;
  elseif ~isnan(after_dist)
      if ~isempty(in_chr)
        closest_idx=in_chr(after_idx(after_min_idx));
        closest_dist=after_dist;
      else
            closest_idx = [];
          closest_dist = 0;
      end
  end
end

    
