function [D focal_regs params] = run_gistic20(base_dir,D,refgene,params)
%RUN_GISTIC20 run focal and broad GISTIC analyses and save results
%
%  [D focal_regs params] = run_gistic20(BASE_DIR,D,REFGENE,PARAMS)
%
%  RUN_GISTIC20 is a top level GISTIC pipeline function that
%  performs a GISTIC focal analysis on copy-number data stored in a D
%  structure by calling run_focal_gistic. It optionally performs an
%  additional broad analysis and optionally saves the results of the
%  analyses as files.
%
%  BASE_DIR is the output directory where GISTIC analyses are written.
%
%  D is a structure containing prepared segmented copy number data.
%
%  REFGENE is the path to a reference gene file, or the contents of a
%    loaded refgene in a struct
%
%  PARAMS is a structure containing optional parameters:
%   ** PARAMS.run_broad_analysis, if set, runs a broad analysis
%    (default cleared)
%   ** PARAMS.write_gene_files, if set, enables the output of
%    broad, focal, and combined broad/focal data by genes.
%   ** PARAMS.gene_collapse_method - string designating how to collapse marker
%    level data to gene level data for the gene files. 'mean' specifies to
%    use the average copy level of the markers (default); 'median' specifies
%    to  use the median marker value; 'extreme' means use the value furthest
%    from zero.
%   ** PARAMS.ext is a string specifying  a penultimate extension
%    added to saved file names to distinguish between runs. It
%    should start with, but not end with, a period. 
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
%    to be considered significant. The default is .25.
%   ** PARAMS.array_list_file names a tab-delimited file with headers whose
%    first column contains the list of arrays on which to run the analysis.
%    If empty, analyzes all samples in D. 
%   ** PARAMS.do_arbitration is a boolean value: if 1 use arbitration
%    algorithm to distribute score between overlapping peaks (the default);
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

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  %% Set optional parameter defaults
  if ~exist('params','var') || isempty(params)
    params = struct();
  end
  % use standard GISTIC2 defaults for parameters that are not defined
  params = gistic2_param_defaults(params);
  
  %% Check required input files...
  
  % base directory
  if ~exist('base_dir','var') || isempty(base_dir)
    throw(MException('snp:gistic:no_base_dir','Must supply base_dir!'));
  else
    base_dir = add_slash_if_needed(base_dir);  
    if ~exist(base_dir,'dir')
      mkdir(base_dir)
    end
  end
  
  % D-struct input
  if ~exist('D','var') || isempty(D) || ~(ischar(D) || isstruct(D))
    throw(MException('snp:gistic:no_D_struct','Must supply valid D file or struct!'));
  end
  if ischar(D) && ~exist(D,'file') 
    throw(MException('snp:gistic:no_D_struct_file','D-struct file ''$s'' does not exist!',D))
  end

  % reference genome file
  if ~exist('refgene','var') || isempty(refgene)
    throw(MException('snp:gistic:no_refgene','Must supply refgene'));
  else
    % load reference genome, set rg and cyto
    refgene = load_refgene(refgene);
    rg = refgene.gene;
    cyto = refgene.cyto;
  end
  
  %% Run focal analysis
  verbose('Running focal analysis...',20);
  [D,focal_regs,params] = run_focal_gistic(base_dir,D,refgene,params);
  params.thresh = [params.t_amp params.t_del params.t_amp params.t_del];
    
  %% Run broad analysis
  if params.run_broad_analysis
    verbose('Running broad analysis...',20);
    [arm_medians,~,~,~,~,~,~,arm_names,~] = gistic_broad_analysis(base_dir,D,refgene,params);
  end
  
  %% Write out gene tables
  pathbase = [base_dir params.fname];
  extxt = [params.ext '.txt'];

  if params.write_gene_files
    
    % Write arm medians
    if exist('arm_medians','var') && ~isempty(arm_medians)
      write_arm_medians([pathbase 'broad_values_by_arm' extxt],D,arm_medians,arm_names);
    end    

    % Check if D is in log space or copy number space
    if ~isfield(D,'islog') 
      warning('Don''t know if D is in log space or copy number space...');
    elseif D.islog
      D.dat = 2.^(D.dat+1)-2;
      D.islog = 0;
    end
    
    % assign markers to genes
    if ~isfield(rg,'snps')
      rg=add_snps_to_rg(rg,D);
    end

    % write out all_data_by_genes
    verbose('Reducing all data to genes',10);
    collapse_and_write_genes([pathbase 'all_data_by_genes' extxt],D,rg,cyto,params);

    % Write out focal_data_by_genes
    verbose('Reducing focal data to genes',10);
    focals = reconstruct_genomes(D.Qs,struct('broad_or_focal', 'focal',...
                                             't_amp', params.t_amp,...
                                             't_del', params.t_del,...
                                             'broad_len_cutoff', params.broad_len_cutoff,...
                                             'column_to_add', 12,...
                                             'use_segarray', params.use_segarray,...
                                             'rows', length(D.pos),...
                                             'cols',length(D.sdesc) ));
    D1 = D;
    D1.dat = focals.amp+focals.aod-focals.del-focals.doa;
    clear focals
    
    collapse_and_write_genes([pathbase 'focal_data_by_genes' extxt],D1,rg,cyto,params);
    
    % Write out broad_data_by_genes
    verbose('Reducing broad data to genes',10);
    broads = reconstruct_genomes(D.Qs,struct('broad_or_focal', 'broad',...
                                           't_amp', params.t_amp,...
                                           't_del', params.t_del,...
                                           'broad_len_cutoff', params.broad_len_cutoff,...
                                           'column_to_add', 12,...
                                           'use_segarray', params.use_segarray, ...
                                           'rows', length(D.pos),...
                                           'cols',length(D.sdesc) ));
    D1.dat = broads.amp+broads.aod-broads.del-broads.doa;
    clear broads
    collapse_and_write_genes([pathbase 'broad_data_by_genes' extxt],D1,rg,cyto,params);
    clear D1
    
    %% Write gene calls
    gene_calls(base_dir,params,refgene);
    verbose('Finished writing gene tables.',10);
  end
  verbose('GISTIC2 run completed!',10);
  
%% subfunction to collapse and write out files
function collapse_and_write_genes(fname,D,rg,cyto,params)
    G = reduce_to_genes(D,rg,'symb', struct('collapse_method',params.gene_collapse_method,'find_snps_type',1));
    write_gene_table(fname,G,cyto,rg);

  
