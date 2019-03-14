function approx_range = interpolate_range_dist_percentile(perm_ads,score_type,conf_level,win_range,interpl_method)
% INTERPOLATE_RANGE_DIST_PERCENTILE find delta-G score as a function of window size
%
%   approx_range = interpolate_range_dist_percentile( ...
%                   perm_ads,score_type,conf_level,win_range,interpl_method)
%
% PERM_ADS a 2xN array of score vectors made from events permuted across the genome
% SCORE_TYPE structure whose 'res' field gives the resolution to use in
% constructing empirical ditributions
% CONF_LEVEL confidence level for the maximum delta-G
%
% TODO - finish documenting this function

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if ~exist('interpl_method','var') || isempty(interpl_method)
    interpl_method = 'linear';
  end
  
  interpl_base = 2;
  max_power = ceil(log2(max(win_range))/log2(interpl_base));
  interpl_vals = interpl_base.^(0:max_power);   
  approx_range = cell(1,length(perm_ads));
  
  for i=1:length(perm_ads)
    disp(i)
    perms = perm_ads{i};
    gt = zeros(length(perms),length(interpl_vals));
    for j=1:length(perms)
      verbose(['Calculating on permutation %d of %d'],30,j,length(perms));
      R = cell(length(interpl_vals),1); % score range for window size
      M = cell(length(interpl_vals),1); % minimum G score for window size
      m = cell(length(interpl_vals),1); % maximum G score for window size
      % initialize for window size 1
      M{1} = perms{j};
      m{1} = perms{j};
      R{1} = M{1}-m{1};
      % create normalized empirical ditribution
      qq = histc(R{1},0:score_type.res:max(R{1}));
      qq = qq/sum(qq);
      gt(j,1) = snp_inv_cdf(qq',score_type,conf_level);
      for k=2:length(interpl_vals)
        %verbose(['Calculating range dist for window size: %d'],30,interpl_vals(k))
        len = interpl_vals(k);
        M{k} = max([M{k-1}(1:(end-len/2)); M{k-1}((len/2+1):end)], ...
                   [],1);
        m{k} = min([m{k-1}(1:(end-len/2)); m{k-1}((len/2+1):end)], ...
                   [],1);
        R{k} = M{k}-m{k};
        % create normalized empirical distribution
        qq = histc(R{k},0:score_type.res:max(R{k}));
        qq = qq/sum(qq); % normalize distribution
        gt(j,k) = snp_inv_cdf(qq',score_type,conf_level);
      end
    end
    % interpolate delta-G for window sizes in between powers of two (log scaled)
    approx_range{i} = interp1(log2(interpl_vals),mean(gt,1),log2(win_range),interpl_method);
  end
  
  
  
    
    
