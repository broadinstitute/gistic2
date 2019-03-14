function [genomes segs] = reconstruct_genomes(Qs,params)
% RECONSTRUCT_GENOMES - create copy number arrays from underlying events
%   [GENOMES SEGS] = reconstruct_genomes(Qs,PARAMS)
%   Qs is a structure arrays containing all SCNA events
%   PARAMS is a structure containing optional parameters
%     PARAMS.broad_len_cutoff defines the threshold below which events
%       are considered focal (versus broad)
%     PARAMS.broad_or_focal is either 'broad' or 'focal' depending on
%       whether the genome should be reconstructed from broad or focal events 
%     PARAMS.t_amp and PARAMS.t_del are amplification and deletion
%       significance thresholds below which events will not be used in genome
%       reconstruction. Default 0.1.
%     PARAMS.column_to_add specifies the column of Qs to use to reconstruct
%       the genome: 4 for score or 12 for event amplitude. The default is 4.
%     PARAMS.rows optionally specifiy the number of samples in the returned
%     genome.
%     PARAMS.cols optionally specifies the number of markers in the
%     returned genome.
%     PARAMS.use_segarray is set to true to use the SegArray compression
%     scheme on the data arrays to same memory
%   
%   The returned GENOMES is a structure of genome position-by-sample copy
%   number arrays with members 'amp', 'del', 'aod', 'doa' reconstructed from
%   the same-nasmed kinds of events in Qs. SEGS is a structure of arrays
%   containing subsets of the input segments corresponding to the Qs,
%   filtered according to the options specified in params.
%
%  See also Qs.

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  params = impose_default_value(params,'broad_or_focal','focal');
  params = impose_default_value(params,'t_amp',0.1);
  params = impose_default_value(params,'t_del',0.1);
  params = impose_default_value(params,'broad_len_cutoff',0.98);
  params = impose_default_value(params,'column_to_add',4);
  params = impose_default_value(params,'use_segarray',true);
  params = impose_default_value(params,'rows',NaN);
  params = impose_default_value(params,'cols',NaN);
  
  fields = intersect({'amp','del','aod','doa'},fieldnames(Qs));

  %% figure out the size of the output genomes CN array

  % set genome length 
  if isnan(params.rows)
    % find maximum segment end (column 3) in all the fields of Qs
    dat = [];
    for i = 1:length(fields)
      dat = [dat;Qs.(fields{i})(:,3)];
    end
    s = max(dat);
  else
    s = params.rows;
  end
  % set sample width
  if isnan(params.cols)
    % find maximum sample (column 5) in all the fields of Qs
    dat = [];
    for i = 1:length(fields)
      dat = [dat;Qs.(fields{i})(:,5)];
      n = max(dat);
    end
  else
    n = params.cols;
  end

  % create empty genome CN template
  if params.use_segarray
    blank = SegArray.zeros(s,n);
  else
    blank = zeros(s,n);
  end
  
  genomes = struct();
  segs = struct();
 
  % loop over genomes being reconstructed
  for j=1:length(fields)
    field = fields{j};
    if strmatch('a',field)
      thresh = params.t_amp;
    else
      thresh = params.t_del;
    end
    
    verbose(['Reconstructing genome: ' field],20);
    
    cur_Q = Qs.(field);
    genomes.(field) = blank;
    
    % select focal or broad segments
    if isequal(params.broad_or_focal,'focal')
      cur_segs = intersect(find(cur_Q(:,8) < params.broad_len_cutoff),find(cur_Q(:,12) ...
                                                        >= thresh));
    else
      cur_segs = intersect(find(cur_Q(:,8) >= params.broad_len_cutoff),find(cur_Q(:,12) ...
                                                        >= thresh));
    end
    segs.(field) = cur_segs;
    
    % accumulate the selected segments into a copy number array
    if params.use_segarray
      % SegArray genome: use addSegments
      genomes.(field) = addSegments(genomes.(field),...
        cur_Q(cur_segs,2),... % start row
        cur_Q(cur_segs,3),... % end row
        cur_Q(cur_segs,5),... % column (sample)
        cur_Q(cur_segs,params.column_to_add)); % values
    else
      % normal array genome: loop over samples
      cur_genome = genomes.(field);
      for i=1:length(cur_segs)
        if mod(i,100) == 0
          verbose([num2str(i) ' of ' num2str(length(cur_segs))],30)
        end
        cur_genome(cur_Q(cur_segs(i),2):cur_Q(cur_segs(i),3),cur_Q(cur_segs(i),5)) = ...
              cur_genome(cur_Q(cur_segs(i),2):cur_Q(cur_segs(i),3),cur_Q(cur_segs(i),5)) + ...
              cur_Q(cur_segs(i),params.column_to_add);
      end
      genomes.(field) = cur_genome;
    end
  end

 
