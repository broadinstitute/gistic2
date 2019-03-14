function [B chrnEnd chrarms] = make_sample_B(D,idx,cyto,chrnEnd,chrarms)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  use_segarray = isa(D.dat,'SegArray');

  if exist('breakpt_ids','var') && exist('sample_ids','var')
    if length(breakpt_ids) ~= length(sample_ids)
      error('Breakpoint and sample ids must have same length!');
    end
  end
  
  if ~exist('chrnEnd','var') || isempty(chrnEnd)
    changeChrn = diff(D.chrn);
    chrnEnd = [find(changeChrn == 1); size(D.chrn,1)];
  end
  
  if ~exist('chrarms','var')
    chrarms = [];
  end
  
  if ~exist('idx','var') || isempty(idx)
    error('Must supply a sample index to make B!');
  end

  % find the CN breakpoints in the data
  if use_segarray
    bpt = find(getbpts(D.dat(:,idx)));
    bpt = bpt(2:end)-1;
  else
    bpt = find(diff(D.dat(:,idx)) ~= 0);
  end
  % add in chromosomal breakpoints
  bpt = union(bpt,chrnEnd);
  % make the "B" array
  samples = repmat(idx,length(bpt),1);
  chrn = full(D.chrn(bpt));
  st = [1; bpt(1:end-1)+1];
  en = bpt(:);
  amp = full(D.dat(bpt,idx));
  
  B = [chrn st en amp samples];
  [B(:,6) chrarms] = normalize_by_arm_length(D,B,cyto,1,2,chrarms);
  
