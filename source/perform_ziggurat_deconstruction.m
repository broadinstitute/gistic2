function [D,Qs] = perform_ziggurat_deconstruction(D,cyto,ziggs,cap)
%PERFORM_ZIGGURAT_DECONSTRUCTION perform the ziggurat algorithm on CN data
%
%   D = perform_ziggurat_deconstruction(D,cyto,ziggs,cap)
%  
%   PERFORM_ZIGGURAT_DECONSTRUCTION adds a 'Qs' copy number event structure 
% to a D structure with input copy number data. The CN values in D.dat are 
% assumed to be in log ratio form unless the 'islog' field is present and set
% to 0 (false).
%
%   CYTO is a structure containing physical chromosome information. 
%   ZIGGS is a structure of optional parameters:
%       ZIGGS.max_segs_per_sample is the segment limit for filtering out
%         noisy samples
%       ZIGGS.seg_count_file specifies an output file reporting segment
%         counts and sample inclusion status
%   CAP is an optional range limit to apply to CN data, a length two vector
%       [min,max]
%
% The Qs structure added to D has four fields for different kinds of events: 
% 'amp' contains pure amplification (over amplification) events, 'del' pure 
% deletion events, 'aod' amplifications over deletions, and 'doa' deletions
% over amplifications. Each field of Qs is an Nx12 array, where N is
% the number of events (rows). The 12 columns of each field hold event
% properties:
%
%   column 1 - chromosome
%   column 2 - start SNP index
%   column 3 - end SNP index
%   column 4 - event amplitude
%   column 5 - sample index (column in D.dat)
%   column 6 - start CN level
%   column 7 - end CN level
%   column 8 - event length (as chromosome arm fraction)
%   column 9 - deconstruction score
%   column 10 - arm score
%   column 12 - event amplitude
% 

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  use_segarray = isa(D.dat,'SegArray');
    
  % Check that required files are present, otherwise error
  varlist1 = {'D','cyto'};
  for idx=1:length(varlist1)
    if ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
      error('Required input %s undefined.',varlist1{idx})
    end
  end
  % process optional arguments
  if ~exist('ziggs','var')
    ziggs = struct;
  end
  if ~isfield(ziggs,'max_segs_per_sample') 
    max_segs_per_sample = 500;
  else
    max_segs_per_sample = ziggs.max_segs_per_sample;
  end

  if ~isfield(ziggs,'seg_count_file') 
    seg_count_file = '';
  else
    seg_count_file = ziggs.seg_count_file;
  end

  % The 'islog' field was a retrofit to the D struct, so the assumption is
  % that if 'islog' is not a D field, the D.dat CN data is log-ratio scaled
  % as it usually is.
  if ~isfield(D,'islog')
     verbose('No ''islog'' field - assuming data are log ratio.',30);
     D.islog = 1;
  end
  
  % cap (range limit) the data
  if exist('cap','var') && ~isempty(cap)
      if isscalar(cap)
          maxcap = cap;
          mincap = -cap;
      else
          % allow "unbalanced" vector cap
          maxcap = cap(1);
          mincap = cap(2);
      end
      if ~D.islog
          maxcap = 2^(1+maxcap)-2;
          mincap = 2^(1+mincap)-2;
      end
      if use_segarray
          D.dat = cap_vals(D.dat,[maxcap mincap]);
      else
          D.dat(D.dat > maxcap) = maxcap;
          D.dat(D.dat < mincap) = mincap;
      end
  end
  
  %% remove noisy samples
  D = remove_noisy_samples(D,max_segs_per_sample,seg_count_file);
  if size(D.dat,2) < 1
    throw(MException('snp:perform_ziggurat_deconstruction:all_data_removed', ...
                    'All samples were removed by noise filtering.'));
  end

  if size(D.dat,2) == 1
    throw(MException('snp:perform_ziggurat_deconstruction:one_sample_left', ...
                    'Only one sample remains after noise filtering.'));
  end
  
  %% Convert to copy number space...
  if D.islog
    verbose('Converting data from log to copy number space...',30);
    D.dat = 2.^(D.dat+1)-2;
    D.islog = 0;
  else
    verbose('Data already in copy number space...',30);
  end

  %% call inner deconstruction function
  verbose('Performing ziggurat deconstruction!',20);    
  [D,QA,QD,QAOD,QDOA] = perform_deconstruction(D,cyto,[],1);
  % massage the data
  QD(:,4) = -1*QD(:,4);
  QDOA(:,4) = -1*QDOA(:,4);
  QA(:,12) = QA(:,4);
  QD(:,12) = QD(:,4);
  QAOD(:,12) = QAOD(:,4);
  QDOA(:,12) = QDOA(:,4);
  % create Qs structure and add it to D
  Qs = struct('amp',QA,'del',QD,'aod',QAOD,'doa',QDOA);
  Qs.header = {'chrn','pos_start (snp)','pos_end (snp)','scna_score',...
               'sample_id','starting_cn_level','ending_cn_level','fract_chr_arm_length',...
               'deconstruction_score','arm_score','EMPTY','amplitude'};
  D.Qs = Qs;
  Qs.sdesc = D.sdesc;
  
