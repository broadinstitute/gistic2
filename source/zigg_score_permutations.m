function [q,p,d,ads,min_val] = zigg_score_permutations(Z,res,cleaned_snps)
%ZIGG_SCORE_PERMUTATIONS permute marker G-scores to get p- and q-values
%
%   [Q,P,D,ADS,MIN_VAL] = ZIGG_SCORE_PERMUTATIONS(Z,RES,CLEANED_SNPS)
%
% in:
%       Z = matrix of amplification or deletion scores
%       RES = resolution used for empirical distributions
%       CLEANED_SNPS = an optional vector of SNP counts for each sample to
%       be removed from the low-scoring end
% out:
%       P = vector of p-values
%       Q = vector of q-values
%       D = sum-of-scores distribution expected by chance
%       ADS = vector of amplification or deletion G-scores
%       MIN_VAL = minimum value in the input score matrix

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


  n = size(Z,2);  % number of samples
  
  if ~exist('cleaned_snps','var') || isempty(cleaned_snps)
    cleaned_snps = zeros(1,n);
  end
  
  % adjust scores so that all are positive
  min_val = min(min(Z,[],1));
  if min_val < 0 
    Z = Z - min_val;
  end
  % compute G-score
  ads = nanmean(Z,2);
  % create array of score histograms, one for each sample
  ha = cell(1,n);
  for i=1:n
    ha{i} = histc(Z(:,i)/n,0:res:(max(Z(:,i))/n+res));
    ha{i} = clean_hist(ha{i},cleaned_snps(i));
    ha{i} = ha{i}/sum(ha{i});
    if mod(i,100)==0
      disp(i);
    end
  end
  
  % calculate expected sum-of-scores null distribution using convolutions
  max_bin = max(1,ceil(max(ads)/res));
  fprintf(1,'Conv:.');
  d = ha{1};
  for i=2:length(ha);  
    d=conv(d,ha{i});
    % move any overflow to edge of distribution window 
    if length(d) > max_bin
      d(max_bin) = sum(d(max_bin:end));
      d = d(1:max_bin);
    end
    fprintf(1,'.');
  end
  fprintf(1,'\n');
  
  % normalize null distribution
  sd = sum(d);
  if abs(sd-1) > 0.1
    error('not a distribution'); % (should only be small errors)
  else
    d = d ./ sd;
  end
  
  % create a reverse cumulative distribution function
  t = flipud(cumsum(flipud(d)));
  % look up the p values from the scores & calculate false discovery rate
  p = t(min(1+floor(ads ./ res),length(t)));
  q = calc_fdr_value(p);
  
