function new_ads = permute_seg_locations(Q,nsamples,num_snps)
% permute_seg_locations -- Permutes segment locations in Q to
% create null scores for computation of range distribution
% 
%
% new_ads = permute_seg_locations(Q,nsamples,nsnps)
%    
%
%
% Q (requird): nseg x 12 matrix containing ziggurat segments for a
%              dataset.  Each row represents a segment to be permuted.
% nsamples(required) = number of samples in dataset
% num_snps(required) = number of snps in dataset
%
% Outputs: new_ads, a 1 x num_snps array containing a permuted score
% distribution calculated by randomly shifting the location of each
% segment in Q across the genome.
%
% ---
% $Id$
% $Date: 2014-01-31 15:30:24 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
    
  shift = floor(num_snps*rand(size(Q,1),1));
  temp_Q = Q;
  temp_Q(:,2:3) = mod(temp_Q(:,2:3)+repmat(shift(:,1),1,2),num_snps)+1;
  cur_ads = zeros(1,num_snps);
  for j=1:size(temp_Q,1)
    if temp_Q(j,2) < temp_Q(j,3)
      cur_ads(temp_Q(j,2):temp_Q(j,3)) = cur_ads(temp_Q(j,2):temp_Q(j, ...
                                                          3))+temp_Q(j,4)/nsamples;
    else
      cur_ads(temp_Q(j,2):num_snps) = cur_ads(temp_Q(j,2):num_snps)+ ...
          temp_Q(j,4)/nsamples;
      cur_ads(1:temp_Q(j,3)) = cur_ads(1:temp_Q(j,3))+temp_Q(j,4)/ ...
          nsamples;
    end
  end
  
  new_ads = cur_ads;
  
  
  
