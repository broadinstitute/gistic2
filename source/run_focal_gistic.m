function [D,regs,params] = run_focal_gistic(base_dir,D,refgene,params)
% RUN_FOCAL_GISTIC -- Perform Focal GISTIC2.0 pipeline
%
% D = run_focal_gistic(BASE_DIR,D,REFGENE,PARAMS)   
%
% Required inputs:
%     BASE_DIR: directory in which to write output for this run
%     D: either file pointing to a saved D_struct, or D_struct itself
%     REFGENE: refgene/cytoband file location, or its contents packed in a
%        struct
% PARAMS is an optional struct containing focal gistic runtime parameters
% Accepted fields to the PARAMS struct include:
%    -- ziggs: structure controlling ziggurat deconstruction. Accepts two fields:
%         1) max_segs_per_sample (default 2500): Remove samples with more than this
%         number of segments
%    -- t_amp (default: 0.1): Amplitude cutoff for amp segments
%    -- t_del (default: 0.1): Amplitude cutoff for del segments
%    -- broad_len_cutoff (default: 0.98): Length cutoff for calling focal
%                                  segments (expressed as fraction of chr arm)
%    -- conf_level (default: 0.75): The confidence level(s) of a peak containing a driver used
%                                  for setting the peak boundaries.If passed a vector of values, 
%                                  multiple peaks for each region will be reported , one for 
%                                  each confidence level.  
%    -- qv_thesh (default = 0.25): Q-value threshold for significance
%    -- res (default = 0.05/# samples): resolution for permutations analysis
%    -- cap (default = 1.5): max absolute log2 copy value
%    -- alpha (default = empty): The alpha scaling parameter to use for scoring amplitude of
%                                segments.  If empty, is calculated from data.
%    -- do_gene_gistic (default = 0): Perform gene gistic on deletions  
%    -- array_list_file (default empty): List of arrays on which to run
%                                        the analysis. If empty, uses all samples in D.
%    -- save_data_files (default = 1): Save intermediate data files
%    -- conserve_disk_space (default = 0): Don't save big data files.
%    -- ext (default = empty): Text extension added to output file names
%    -- do_arbitration (default = 1): Run Arbitrated peel-off (1) or
%                                     regular greedy peel-off (0)
%    -- arm_peeloff (default = 0):   Peel off just segments (0) or segments
%                                     along with other events on the same
%                                     arm for a given sample
%    -- peak_types (default = 'robust'): Type of wide peaks to
%                                        calculate. Accepted values are
%                                        'robust' and 'loo' (leave one
%                                        out).  If both given, whichever type is listed first
%                                        is used for gistic output files.             

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
  %% Set defaults
  
  if ~exist('params','var') || isempty(params)
    params = struct();
  end
  % use standard GISTIC2 defaults for parameters that are not defined
  params = gistic2_param_defaults(params);
    
  %% Check required variables
  if ~exist('base_dir','var') || isempty(base_dir)
    throw(MException('snp:gistic:no_base_dir','Must supply base_dir!'));
  else
    base_dir = add_slash_if_needed(base_dir);  
    if ~exist(base_dir,'dir')
      mkdir(base_dir)
    end
  end

  % handy string abbreviations
  pathbase = [base_dir params.fname]; % path + base name
  ext = params.ext;                   % user extension

  % D-struct input
  if ~exist('D','var') || isempty(D) || ~(ischar(D) || isstruct(D))
    throw(MException('snp:gistic:no_D_struct','Must supply valid D file or struct!'));
  end
  if ischar(D) 
    % D is path to saved D-struct
    if exist(D,'file')
        D = load_D(D);
    else
        throw(MException('snp:gistic:no_D_struct_file','D-struct file ''$s'' does not exist!',D))
    end
  else
    verbose('Using supplied D...',20);
  end
  
  params.use_segarray = isa(D.dat,'SegArray');

  % load reference genome
  if ~exist('refgene','var') || isempty(refgene)
    throw(MException('snp:gistic:no_refgene','Must supply refgene'));
  else
    refgene = load_refgene(refgene);
    rg = refgene.gene;
    cyto = refgene.cyto;
  end
    
  %% Subselect samples for analysis 
  if ~isempty(params.array_list_file)
    if ~iscell(params.array_list_file)
      AL = read_array_list_file(params.array_list_file);
      array_list = {AL.array};
    else
      array_list = params.array_list_file;
    end
    verbose(['Subselecting ' num2str(length(array_list)) ' samples in data.'],30);
    [~,~,use_arrays]=intersect(array_list,D.sdesc);
    D = reorder_D_cols(D,use_arrays);
    write_array_list([pathbase 'arraylistfile' ext '.txt'],D);
  end
  
  % adjust resolution by sample size
  Nsamples = size(D.dat,2);
  params.res = .05/Nsamples;
  score_type.res = params.res;
  score_type.cap = params.cap;
  
  % limit amplification score based number of bins that may be required for sample
  % (NOTE: overall copy number used as limiting proxy for actual score)
  MAX_SAMPLE_BINS = 1e6;
  limitZ = MAX_SAMPLE_BINS*Nsamples*params.res*max(params.alpha);
  if ~isfield(D,'islog') || D.islog
      limitZ = log2(limitZ+2)-1;
  end
  too_high = max(D.dat) > limitZ;
  if any(too_high)
      warning('gistic:memory:score_too_high',...
              'Excessive score truncated for %d samples.',sum(too_high));
      D.dat(D.dat > limitZ) = limitZ;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% BEGIN GISTIC 2.0 pipeline
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  verbose('Running focal GISTIC version %s',10,gistic_version);
  verbosedisp(params,10);
  
  
  %% STEP 1: Ziggurat deconstruction

  % if needed, perform ziggurat deconstruction
  if ~isfield(D,'Qs')
    params.ziggs.seg_count_file = [pathbase,'sample_seg_counts',ext,'.txt'];
    D = perform_ziggurat_deconstruction(D,cyto,params.ziggs,params.cap);
  end
  
  % save segmented data and deconstructed events
  if params.save_data_files && ~params.conserve_disk_space
    save_D([pathbase 'D.cap' num2str(params.cap) ext '.mat'],D,'-v7.3');
  end
  
  [nmarkers,nsamples] = size(D.dat);

  %% STEP 2: Reconstruct focal genome

  verbose('Reconstructing focal genomes...',20);
    
  [focals,focal_segs] = reconstruct_genomes(D.Qs,struct('broad_or_focal','focal',...
                                                    'broad_len_cutoff',params.broad_len_cutoff,...
                                                    't_amp',params.t_amp,...
                                                    't_del',params.t_del,...
                                                    'column_to_add',12,...
                                                    'use_segarray',params.use_segarray,...
                                                    'rows',nmarkers,...
                                                    'cols',nsamples) );
  % write binary file of focal segments
  if params.save_data_files && ~params.conserve_disk_space
    focal_out = [pathbase 'focal_dat.' num2str(params.broad_len_cutoff) ext '.mat'];
    verbose('Saving results to: %s',20,focal_out);
    save(focal_out,'focals','focal_segs','-v7.3');
  end
  
  % write text file of focal segments
  if ~params.genepattern
    Z = D;
    Z.dat = focals.amp-focals.del;
    write_seg_file([pathbase 'focal_input' ext '.seg.txt'],Z)
    clear Z
  end
  
  %% STEP 3: Calculate genome-wide scores

  [scores,D,params] = score_genome(D,focals,params);
  clear focals
  
  if params.save_data_files && ~params.conserve_disk_space
    score_out = [pathbase 'scores.' num2str(params.broad_len_cutoff) ext '.mat'];
    verbose('Saving results to: %s',20,score_out);
    save(score_out,'scores','-v7.3');
  end
  
  %% STEP 4: Perform background permutations 
  
  [q,p,d,ads,score_thresh,gg_rg,gg_ds] = compute_stats(D,scores,rg,score_type,...
                                                    params.alpha,params.qv_thresh,params.do_gene_gistic);
  % compute_stats() can veto gene-gistic by returning empty gg_rg
  if isempty(gg_rg)
      params.do_gene_gistic = false;
  end
  
  if params.save_data_files && ~params.genepattern
    verbose('Saving stats...',20);
    t_amp = params.t_amp; t_del = params.t_del;
    alpha = params.alpha;
    save([pathbase 'orig_stats' ext '.mat'],...
                'ads','d','p','q','score_thresh','alpha','t_amp','t_del');
    
    if params.do_gene_gistic
      save([pathbase 'gene_stats' ext '.mat'],...
                'gg_rg','gg_ds');
    end
  end
  
  %% STEP 5: Run arbitrated peel-off
  
  verbose('Running arbitrated peel-off...',20);
  Z.dat{1} = scores.amp-min(min(scores.amp)); Z.dat{2} = scores.del-min(min(scores.del));
  Z.pos = D.pos;
  Z.chrn = D.chrn;
  % add arm information to scores if needed
  if params.arm_peeloff
    Z = add_cyto(Z,cyto);
  end
  clear scores
  
  % focQs = focal pure amp and del segments
  focQs = D.Qs;
  focQs = rmfield(focQs,{'aod','doa'});
  focQs.amp = focQs.amp(focal_segs.amp,:);
  focQs.del = focQs.del(focal_segs.del,:);
  clear focal_segs
  
  [regs,seg_assignment] = identify_peaks_by_arbitration(Z, focQs, ads, q, score_thresh, ...
                            gg_rg,gg_ds,score_type.res*params.alpha(2),...
                            params.do_arbitration,params.arm_peeloff);
  
  if params.save_data_files
    verbose('Saving peak regions...',20);
    save([pathbase 'peak_regs' ext '.mat'],'regs','seg_assignment');
  end
  
   
  %% STEP 6: ROBUST
  
  % Compute permuted score distributions
    
  verbose('Computing permuted score distributions...',20);
  nperm = 25;
  
  if params.do_gene_gistic
    clean_struct = struct('clean_up',1,'seg_assignment',seg_assignment, ...
                          'score_type',score_type,'d',{d{1} gg_ds{1}},'alpha',params.alpha);
  else
    clean_struct = struct('clean_up',1,'seg_assignment',seg_assignment, ...
                          'score_type',score_type,'d',d,'alpha',params.alpha);
  end

  perm_ads = permute_segment_locations(focQs,nperm,nsamples,nmarkers,clean_struct);
  
  if params.save_data_files
    verbose('Saving permuted score files...',20);
    save([pathbase 'perm_ads' ext '.mat'],'perm_ads');
  end
  
  % Compute wide peaks
  verbose('Computing wide peak regions...',20);
  
  range_struct = struct('method','approximate','interpl_method','linear');
  
  params.conf_level = unique(params.conf_level);
  ncl = length(params.conf_level);
  wregs = cell(1,ncl);
  
  % loop over confidence values to determine boundaries
  for l=1:ncl
    verbose('Running Regbounder at confidence level: %d%%',30,100*params.conf_level(l))
    wregs{l} = calculate_wide_peaks_for_regs(Z,focQs,regs,perm_ads,params.conf_level(l),...
                                               score_thresh,score_type,range_struct,...
                                               gg_rg,params.peak_types);
   end
  
  % save all wide peaks in matlab format
  if params.save_data_files
    verbose('Saving final_regions...',20);
    save([pathbase 'wide_peak_regs' ext '.mat'],'wregs');
  end
  
  % End ROBUST
  
  %% STEP 7: Save output files
  
  verbose('Writing output files...',20);
  
  % setup partial_hits for amps/dels, which if true for the SCNA type
  % will include genes which straddle the peak boundary in the peak
    
  %! 4/21/2014 now params.partial_hits w/defaults
  %!if params.do_gene_gistic
  %!  partial_hits = [1 0];
  %!else
  %!  partial_hits = [1 1];
  %!end

  ts = [params.t_amp params.t_del];
  
  % loop over confidence values to output files that depend on them
  for l=1:ncl
    cur_regs = wregs{l};
    cext = ['.conf_' num2str(100*params.conf_level(l)) ext];
    
    % Write GISTIC output files
    write_gistic_outfiles(D,base_dir,cext,cyto,cur_regs,ts,...
                          q,p,q,ads,rg,1,score_thresh,0,params.genepattern,...
                          params.partial_hits,0,params.fname);
    % write IGV/UCSC browser file for the regions
    write_bed_file([pathbase 'regions_track' cext '.bed'],...
                   ['GISTIC' cext],'0,0,255',cur_regs,D,cyto);
  end
  
  %% create gene gistic marker scores/q values
  if params.do_gene_gistic
    q{2} = interpolate_gene_scores(D,gg_rg,'q',@min);
    ads{2} = interpolate_gene_scores(D,gg_rg,'gene_scores',@max);
  end
    
  %% save tab-delimited file of scores
  verbose('Writing q value, p value and amp/del score file',30);
  if length(D.pos)==length(q{1})
    write_score_file([pathbase 'scores' ext '.gistic'],D,p,q,ads,ts);
  else
    warning('Length of struct (%d) does not match length of scores (%d)!',...
                size(D.pos,2),length(q{1}));
  end

  %% Create plots of score/q-values and raw data
  
  % generate the plots
  %! TODO!!! need to ensure no real dependency on confidence-sensitive fields of cur_regs)
  all_lesions_file = [base_dir 'all_lesions_file' cext '.txt'];
  gistic_plots(base_dir,ext,D,q,ads,cur_regs,cyto,all_lesions_file,[],[],[],[],[],params.qv_thresh,...
                          [],[],[],0,params.genepattern,1,1,[],[],params.fname);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% End of focal GISTIC 2.0 pipeline
  verbose('Focal GISTIC completed without error',10);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
