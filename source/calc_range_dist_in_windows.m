function range_dist = calc_range_dist_in_windows(bg_dists,win_size,res)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  ndists = size(bg_dists,2);
  nsnps = size(bg_dists{1},2);
  
  tot_range = [];
  for j=1:ndists
    verbose('Calculating on permutation %d  of %d',50,[j ...
                        ndists]);
    cur_dist = bg_dists{j}; 
    cur_range = zeros(1,nsnps-win_size);
    max_diff = max(cur_dist)-min(cur_dist);
    cur_range = zeros(1,nsnps-win_size+1);
    [cur_max cur_max_idx] = max(cur_dist(1:win_size));
    [cur_min cur_min_idx] = min(cur_dist(1:win_size));
    cur_range(1) = cur_max-cur_min;
    for i=2:nsnps-win_size+1
      if i>cur_max_idx %% cur_max was at start of last window
        [cur_max cur_max_idx] = max(cur_dist(i:i+win_size-1));
        cur_max_idx = cur_max_idx + i-1;
      end
      if i>cur_min_idx %% cur_min was at start of last window
        [cur_min cur_min_idx] = min(cur_dist(i:i+win_size-1));
        cur_min_idx = cur_min_idx + i-1;
      end
      if cur_dist(i+win_size-1) > cur_max
        cur_max = cur_dist(i+win_size-1);
        cur_max_idx = i+win_size-1;
      elseif cur_dist(i+win_size-1) < cur_min
        cur_min = cur_dist(i+win_size-1);
        cur_min_idx = i+win_size-1;
      end
      cur_range(i) = cur_max-cur_min;
      
    end
    tot_range = [tot_range cur_range];
  end
  range_dist = histc(tot_range,0:res:max(tot_range));
  range_dist = range_dist/sum(range_dist);
  
