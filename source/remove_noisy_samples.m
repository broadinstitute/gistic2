function D = remove_noisy_samples(D,max_segs_per_sample,seg_count_file)
  % Removes samples with greater than max_segs_per_sample number of
  % segments (default 500)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
  if ~exist('max_segs_per_sample','var') || isempty('max_segs_per_sample')
    max_segs_per_sample = 500;
  end
  % determine number of segments for eac sample
  if isa(D.dat,'SegArray')
        num_bpts = getbpt_counts(D.dat);
  else % uncompressed
      [~,sample_ids] = find(diff(D.dat) ~=0);
      num_bpts = zeros(1,size(D.dat,2));
      for i=1:size(D.dat,2)
        num_bpts(i) = length(find(sample_ids==i));
      end
  end
  % determine samples with acceptable segment counts
  keepers = num_bpts <= max_segs_per_sample;
  % optionally write output file with segment counts
  if exist('seg_count_file','var') && ~isempty(seg_count_file)
      ny = {'no';'yes'};
      segcounts = struct('sample',D.sdesc(:),'segment_count',num2cell(num_bpts(:)),'included',ny(keepers(:)+1));
      savestruct(segcounts,seg_count_file);
  end
  % remove samples with too many segments
  Nremoved = size(D.dat,2) - sum(keepers);
  D = reorder_D_cols(D,keepers);
  verbose('Removed %d samples with more than %d segments',20,Nremoved,max_segs_per_sample);
  
    
