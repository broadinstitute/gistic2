function params = gistic2_param_defaults(params)
%GISTIC2_PARAM_DEFAULTS set standard optional parameter defaults for GISTIC 2.x
%
%       PARAMS = gistic2_param_defaults(PARAMS)
%
%  Processes the structure PARAMS, whose fields are optional GISTIC2
%  parameters. If a given parameter does not have a PARAMS field, the field
%  is added with a default value.
%  
%   ** PARAMS.cnv_file specifies a file containing the genomic locations of
%    germline copy number variants which should be excluded from the analysis.
%   ** PARAMS.join_segment_size specifies the minimum size for a segment:
%    segments smnaller than this value are combined to their largest
%    neighboring segment - the copy number of the resulting segment is a
%    weighted average. The default value is 8.
%   ** PARAMS.remove_X is a  boolean that removes chromosomes numbered
%    above 22 from the analysis. The default is 0.
%   ** use_segarray is a boolean controlling memory performance: if 1, then
%    data is stored in a compressed format that takes more time to access.
%    The default is 1 (as of release 2.0.23).
%   ** PARAMS.run_broad_analysis, if set, runs a broad analysis
%    (default cleared)
%   ** PARAMS.write_gene_files, if set, enables the output of
%    broad, focal, and combined broad/focal data by genes.
%   ** PARAMS.gene_collapse_method - string designating how to collapse marker
%    level data to gene level data for the gene files. 'mean' specifies to
%    use the average copy level of the markers (default); 'median' specifies
%    to  use the median marker value; 'extreme' means use the value furthest
%    from zero, 'max' or 'min' specifies use maximum or minumum respectively.
%   ** PARAMS.ext is a string specifying the penultimate extension
%    added to saved file names to distinguish between runs. It
%    should start with, but not end with, a period. 
%   ** PARAMS.fname is a string specifying the optional base name for the 
%   ** PARAMS.do_gene_gistic, if set, does gene_gistic for the
%    deletion analysis (default cleared)
%   ** PARAMS.t_amp and PARAMS.t_del are log2 copy number values
%    used to filter amplifications and deletions, and are also
%    passed to run_focal_gistic and gistic_broad_analysis.
%   ** PARAMS.broad_len_cutoff is the length cutoff used to
%    distinguish between broad and focal events. 
%   ** PARAMS.conf_level specifies the confidence levels of a peak
%    containing a driver used for setting the peak boundaries. If
%    passed a vector of values, GISTIC2 will report multiple peaks
%    for each region, one for each confidence level.  The default
%    confidence value is .75 (e.g. 75% confidence).
%   ** PARAMS.alpha specifies the exponential [amplification,deletion]
%    factors used to model the frequency dependence on event amplitude. The
%    default values are [2.5145 2.1653].
%   ** PARAMS.cap is a limit that is imposed on the event amplitudes
%    (+/- log ratio units). The default value is 1.5.
%   ** PARAMS.res specifies the resolution used to construct the empirical
%    frequency distributions used to calculate the background frequency.
%    Smaller values are more accurate, but result in longer computation
%    times. The default value is 0.05.
%   ** PARAMS.qv_thresh specifies the minimum q-value that a peak must have
%    to be considered significant. The default is 0.25.
%   ** PARAMS.array_list_file names a tab-delimited file with headers whose
%    first column contains the list of arrays on which to run the analysis.
%    If empty, analyzes all samples in D. 
%   ** PARAMS.do_arbitration is a boolean value: if 1 use arbitration
%    algorithm to distribute score between overlapping peaks (the default);
%    if 0 use greedy algorithm to give overlapping scores to the biggest peak.
%   ** PARAMS.arm_peeloff peel off just segments (0, the default) or segments
%    along with other events on the same arm for a given sample (1)
%    if 0 use greedy algorithm to give overlapping scores to the biggest peak.
%   ** PARAMS.peak_types specifies the type of wide peaks to calculate: 
%    'robust' for robust algorithm, 'loo' for leave-one-out algorithm.
%   ** PARAMS.use_two_sided if true, generate 2D quadrant figure for GISTIC
%    broad analysis. The default is false.
%   ** PARAMS.save_data_files is a boolean which controls the saving of
%    matlab data files to persistant storage. The default is 1.
%   ** PARAMS.conserve_disk_space is a boolean value which further
%    controls the data saving behavior: if 1, than the raw D structure,
%    focal segments and marker scores are not saved to persistant storage.
%   ** PARAMS.ziggs - is a structure containing ziggurat analysis parameters:
%      -- PARAMS.ziggs.max_segs_per_sample is the maximum number of
%        segments a sample can have before it is filtered out of the
%        analysis.
%      -- PARAMS.ziggs.seg_count_file specifies an output file reporting 
%         segment counts and sample inclusion status
%   ** PARAMS.partial_hits - 2-element logical vector (for amp/del) set to
%      include genes that overlap the peak boundary in gene list for the peak
%   ** PARAMS.sample_center - how sample data should be centered, either 
%      'median', 'mean' or 'none'

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


  if ~exist('params','var') || isempty(params)
    params = struct();
  end
  % supply defaults for parameters that are not defined
  params = impose_default_value(params,'t_amp',0.1);
  params = impose_default_value(params,'t_del',0.1);
  params = impose_default_value(params,'broad_len_cutoff',0.98);
  params = impose_default_value(params,'conf_level',0.75);
  params = impose_default_value(params,'qv_thresh',0.25);
  params = impose_default_value(params,'res',0.05);
  params = impose_default_value(params,'cap',1.5);
  params = impose_default_value(params,'alpha',[2.5145 2.1653]);
  params = impose_default_value(params,'do_gene_gistic',0);
  params = impose_default_value(params,'array_list_file','');
  params = impose_default_value(params,'save_data_files',1);
  params = impose_default_value(params,'conserve_disk_space',0);
  params = impose_default_value(params,'fname','');
  params = impose_default_value(params,'ext','');
  params = impose_default_value(params,'use_segarray',1);
  params = impose_default_value(params,'run_broad_analysis',0);
  params = impose_default_value(params,'write_gene_files',0);
  params = impose_default_value(params,'gene_collapse_method','mean');
  params = impose_default_value(params,'use_two_sided',0);
  params = impose_default_value(params,'do_arbitration',1);
  params = impose_default_value(params,'arm_peeloff',0);
  params = impose_default_value(params,'peak_types','');
  params = impose_default_value(params,'genepattern',0);
  params = impose_default_value(params,'sample_center','median');
  % default for partial_hits depends on gene_gistic setting
  if params.do_gene_gistic
    params = impose_default_value(params,'partial_hits',[1,0]);
  else
    params = impose_default_value(params,'partial_hits',[1,1]);
  end
  % fix output base name to have an ending '.'
  if ~isempty(params.fname) && params.fname(end) ~= '.'
      params.fname = [params.fname '.'];
  end
   % fix extension to have a leading '.'
  if ~isempty(params.ext) && params.ext(1) ~= '.'
      params.ext = ['.' params.ext];
  end
  % ziggurat parameter structure
  params = impose_default_value(params,'ziggs',struct);
  params.ziggs = impose_default_value(params.ziggs,'max_segs_per_sample',2500);
  %!params.ziggs = impose_default_value(params.ziggs,'seg_count_file',...
  %!                  [params.fname,'sample_seg_counts',params.ext,'.txt']);
  
  
 