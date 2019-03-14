function snp_segments = find_snp_segments(Q,snp_idx,use_arm)
% find_snp_segments - finds all segments in ziggurat matrix Q covering a
% given snp or range of snps

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  
% snp_segments = find_snp_segments(Q,snp_idx,use_arm)

% Q (required): nsegs x 12 matrix containing ziggurat segments
%               Q(:,2) contains the start snp and Q(:,3) contains
%               the end snp of a segment
% snp_idx (required): snp or set of snps used to identify which
%                     segments to be peeled-off.  Any segment which
%                     covers one or more of the snps in snp_idx will be 
%                     returned.
% use_arm (optional): returned snp_segments will include any segments on
%                     the same arm and same same sample as segments
%                     covering the snps. Arm identifiers are in Q(:,10).
%                     default false
%
% Outputs: snp_segments = n x 1 vector containing indices of snp_segments
% in Q covering the range snp_idx

  if ~exist('use_arm','var') || isempty(use_arm)
      use_arm = false;
  end
  
  % find event segments that overlap the snps
  eventmap = false(size(Q,1),1);
  for j=1:length(snp_idx)
      eventmap = eventmap | (Q(:,2) <= snp_idx(j) & ...
                             Q(:,3) >= snp_idx(j) );
  end
  snp_segments = find(eventmap);
  
  % arm-level peel-off option: add events on same chromosome arm/sample
  if use_arm
      for i=1:length(snp_segments)
          eventmap = eventmap | (Q(:,5) == Q(snp_segments(i),5) & ...
                                Q(:,10) == Q(snp_segments(i),10) );
      end
      snp_segments = find(eventmap);
  end
  
%! original code was heap unfriendly
%{
  if length(snp_idx) > 1
    snp_segments = [];
    for j=1:length(snp_idx);
      snp_segments = [snp_segments;find(Q(:,2) <= snp_idx(j) & Q(:,3) >= ...
                                        snp_idx(j))];
    end
    snp_segments = unique(snp_segments);
  else
    snp_segments = find(Q(:,2)<=snp_idx & Q(:,3)>=snp_idx);
  end
%}
 
    
   
