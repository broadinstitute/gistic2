function [M,Qnew,Qrm] = ziggurat_peeloff(M,Q,idx,snp_or_segs,arm_peeloff)
% ziggurat_peeloff -- Peels off ziggurat segments in Q from a data_matrix
% M that cover a given snp or range of snps
% 
%   [M,Qnew,Qrm] = ziggurat_peeloff(M,Q,snp_idx) 
%    
%
% M (requird): nsnps x nsamples matrix containing the current data matrix
% Q (requird): nsegs x 12 matrix containing ziggurat segments generating M
% idx(required) = set of indices used to identify segments to be peeled-off.
%                 Depending on value of snp_or_segs flag, either
%                 represents the snps to be peeled away or the segments
%                 in Q to be removed 
% snp_or_segs = 1 for snps (default); = 2 for segments): specifies
%                 whether the indices in idx represent snps or segments. 
% arm_peeloff (optional, default = false): if true, segments from the same
%                 sample and chromosome arm are included among those peeled
%                 off
%
% Outputs: M = new_data matrix which results from peeling off segments in
%             Q covering snp_idx
%         Qnew = new 'Q' matrix containing those ziggurat segments remaining after ziggurat peel-off
%         Qrm = 'Q' matrix containing ziggurat segments covering snp_idx

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% ---
% $Id$
% $Date: 2014-01-31 15:31:40 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$  
    
  use_segarray = isa(M,'SegArray');
    
  if ~exist('snp_or_segs','var') || isempty(snp_or_segs)
    snp_or_segs = 1;
  end
   if ~exist('arm_peeloff','var') || isempty(arm_peeloff)
      arm_peeloff = false;
  end

  if snp_or_segs == 1
    snp_segments = find_snp_segments(Q,idx,arm_peeloff);
  else
    snp_segments = idx;
  end
  
  Qrm = Q(snp_segments,:);
  
  verbose('Removing %d segments from data',20,length(snp_segments));
  if use_segarray
      M = addSegments(M,Q(snp_segments,2),Q(snp_segments,3),...
                        Q(snp_segments,5),-Q(snp_segments,4));
  else
      for i=1:length(snp_segments)
        if mod(i,100) == 0
          verbose('Segment %d of %d',30,[i length(snp_segments)]);
        end
        M(Q(snp_segments(i),2):Q(snp_segments(i),3),Q(snp_segments(i),5)) = ...
              M(Q(snp_segments(i),2):Q(snp_segments(i),3),Q(snp_segments(i),5)) ...
              - Q(snp_segments(i),4);
      end
  end
  
  rem_segments = setdiff(1:size(Q,1),snp_segments);
  Qnew = Q(rem_segments,:);
