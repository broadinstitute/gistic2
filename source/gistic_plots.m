function gistic_plots(base_dir,ext,D,q,ads,regs,cyto,all_lesions_file,qv_scale,...
        aLabels,dLabels,aTopLabels,dTopLabels,qt,ull,write_pdf,no_stripes,add_broad_bars,...
        genepattern,plot_raw_data,make_png,flip_dels,gene_name_side,fname)
%Creates GISTIC figures for amplification & deletion using CL21, regs, stats, and cyto
%

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% call from run_focal_gistic:
% gistic_plots(base_dir,cur_ext,D,q,ads,cur_regs,cyto,[],[],[],[],[],[],params.qv_thresh,[],[],[],0,1,1,1);
% default:
% all_lesions_file,qv_scale,aLabels,dLabels,aTopLabels,dTopLabels,
% ull,wp,no_stripes,add_broad_bars, gene_name_side
%
% qt = params.qv_thresh
% add_broad_bars = 0
% plot_raw_data= 1
% make_png = 1
% flip_dels = 1
  
if length(D.pos)~=length(q{1})
  error('Structure and GISTIC output file lengths do not match');
end

if ~exist('regs','var') || isempty(regs)
    regs = {[],[]};
end

if ~exist('flip_dels','var') || isempty(flip_dels)
  flip_dels=0;
end

if ~exist('gene_name_side','var') || isempty(gene_name_side)
  gene_name_side=[1 1-2*(flip_dels>0)]; % 1 is right, -1 is left
end

% q-value significance threshold
if ~exist('qt','var')||isempty(qt)
  qt=0.25;
end

% use loglog scale 1=log10(-log10(q)+1)
if ~exist('ull','var')||isempty(ull)
  ull=1;
end

if ~exist('make_png','var') || isempty(make_png)
  make_png = 0;
end

if ~exist('write_pdf','var') || isempty(write_pdf)
  write_pdf = true;
end

if ~exist('qv_scale','var')
  qv_scale=[];
end

if ~exist('no_stripes','var')
  no_stripes=[];
end

if ~exist('add_broad_bars','var')
  add_broad_bars=[];
end

if ~exist('aLabels','var')
  aLabels = [];
end

if ~exist('dLabels','var')
  dLabels = [];
end

if ~exist('aTopLabels','var')
  aTopLabels = [];
end

if ~exist('dTopLabels','var')
  dTopLabels = [];
end

if ~exist('genepattern','var')
  genepattern = 0;
end

if ~exist('plot_raw_data','var')
  plot_raw_data = 1;
end

if ~exist('fname','var')
  fname = '';
end

% clear any displayed figures
close all

szx = [];

D = shrink_D(D);
if ~isfield(D,'cyton')
    D = add_cyto(D,cyto);
end


%% display raw data
if plot_raw_data
    close('all');
    display_D_heatmap(D);
    % save figure to disk
    
    verbose('About to save Raw_copy_number figure file',30);
    pathname = [base_dir fname 'raw_copy_number' ext];
    
    try
        file_formats = [repmat({'pdf'},write_pdf),...
                        repmat({'png'},write_pdf & make_png),...
                        repmat({'fig'},~genepattern)];
        save_current_figure(pathname,file_formats);
%{
        if ~genepattern
          saveas(gcf,[pathname '.fig'],'fig');
        end
        if write_pdf
          print('-painters','-dpdf',[pathname '.pdf']);
          if make_png
            saveas(gcf,[pathname '.png'], 'png');
          end
        end
%}
    catch me
        verbose('Problems saving raw copy number figures: %s',10,me.message);
    end
end

%% construct reg_names from regions using cytoband

% amplifications
if ~isempty(regs) && size(regs{1},2)>0
    [amp_Peaks_Sorted, amp_IDX] = sort([regs{1}.peak]);
    reg_names{1} = {cyto(D.cyton(amp_Peaks_Sorted)).name};
else
    reg_names{1}=[];
end

% deletions
if ~isempty(regs) && size(regs{2},2)>0
    [del_Peaks_Sorted, del_IDX] = sort([regs{2}.peak]);
    reg_names{2} = {cyto(D.cyton(del_Peaks_Sorted)).name};
else
    reg_names{2}=[];
end

%% get gene_names 

if exist('all_lesions_file','var') && ~isempty(all_lesions_file) && exist(all_lesions_file,'file') 
    gene_names = get_all_lesions_genes(all_lesions_file,regs,reg_names); 
else  %%% gene names are not supplied
    gene_names=reg_names;
end % if we have an all_lesions file

%% create ticks and labels

labels = {aLabels,dLabels};
topLabels = {aTopLabels,dTopLabels};
for k=1:2
    opts = struct;
    if isempty(qv_scale)
        opts.qv_scale = [];
    else
        opts.qv_scale = qv_scale(k);
    end
    opts.qvLabels = labels{k};
    opts.gsLabels = topLabels{k};
    opts.qv_thresh = qt;
    opts.use_loglog = ull;
    [qv_ticks(k),score_ticks(k)] = qplot_scaleticks(q{k},ads{k},opts);
end


%% create regions sorted by position
regs_sort=regs;
for k=1:2 % amps, then dels
  if ~isempty(regs{k})
    pos=cat(1,regs{k}.peak);
    if isfield(regs{k},'added_manually')
       idx=find(cat(1,regs{k}.added_manually)==0);
       pos=pos(idx);
    else
        idx=1:length(pos);
    end
    [~,sposi]=sort(pos);
    regs_sort{k}(idx)=regs_sort{k}(idx(sposi));
  end
end
clear regs

close all

% remove gene names and sorted regions if gene_names begins with '---' why?
del_reg={[],[]};
for k=1:2
  for i=1:length(gene_names{k})
    if strmatch(gene_names{k}{i},'---')
      del_reg{k}=[del_reg{k} i];
    end
  end
  regs_sort{k}(del_reg{k})=[];
  gene_names{k}(del_reg{k})=[];
end

plot_snp_score( base_dir, ext, D, q, ads,qt,ull,write_pdf,'vert',cyto,qv_scale,qv_ticks, ...
                    score_ticks,szx,[1 1],no_stripes,regs_sort,...
                    add_broad_bars,gene_names,genepattern,make_png,flip_dels,gene_name_side,fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SUBFUNCTION - set up scale ticks and threshold
function [qv_ticks,score_ticks] = qplot_scaleticks(q,ads,opts)
%QPLOT_SCALE choose sensible axes for amplification or deletion plot
%
% opts.qv_scale - minimum (left) value for qv_scale
% opts.qvLabels
% opts.gsLabels
% opts.qv_thresh - significance threshold for q-value
%
% optionasl parameter defaults
if ~exist('opts','var') || ~isstruct(opts)
    opts = struct;
end
opts = impose_default_value(opts,'qv_scale',[]);
opts = impose_default_value(opts,'qvLabels',[]);
opts = impose_default_value(opts,'gsLabels',[]);
opts = impose_default_value(opts,'qv_thresh',0.25);
opts = impose_default_value(opts,'use_loglog',3);

% calculate score at threshold
q = q + 1e-323;
sigq = find(q <= opts.qv_thresh);
nsq = find(q > opts.qv_thresh);
if isempty(sigq)
    val = max(ads)+0.1;
else
    % interpolate score value for threshold
    [v0 i0] = max(ads(nsq));
    [v1 i1] = min(ads(sigq));
    lq0 = log10(q(nsq(i0)));
    lq1 = log10(q(sigq(i1)));
    val = v0+(v1-v0)*(log10(opts.qv_thresh)-lq0)/(lq1-lq0);
end

% set scale
if isempty(opts.qv_scale)
    minQ = 0.5*min(q);
else
    minQ = opts.qv_scale;
end

%% create q-value ticks
if isempty(opts.qvLabels)
    if opts.use_loglog == 1
        qv_ticks = loglog_q_ticks(minQ,6,opts.qv_thresh);
    else
        % calculate from minimum Q and significance level
        qvLabel(7) = log10(minQ);
        qvLabel(1) = log10(opts.qv_thresh);
        for k = 2:6
            qvLabel(k) = qvLabel(1) + (qvLabel(7)-qvLabel(1))/(2^(7-k));
        end
        scaleA = 10.^(qvLabel);
        textScaleA = { num2str(scaleA(1)),...
                       ['10^{' num2str(fixdig(qvLabel(2),2)) '}'], ...
                       ['10^{' num2str(fixdig(qvLabel(3),2)) '}'], ...
                       ['10^{' num2str(fixdig(qvLabel(4),2)) '}'], ...
                       ['10^{' num2str(fixdig(qvLabel(5),2)) '}'], ...
                       ['10^{' num2str(fixdig(qvLabel(6),2)) '}'], ...
                       ['10^{' num2str(fixdig(qvLabel(7),2)) '}']  };

        qv_ticks=struct('vals',{scaleA},'txt', {textScaleA});
    end
else
    % q-value ticks provided to the function
    [A1, B1] = get_labels(opts.qvLabels);
    qv_ticks=struct('vals',{A1},'txt', {B1});    
end

%% create score ticks
if isempty(opts.gsLabels)
    % calculate score ticks
    score_ticks = struct('vals',{[val(1) 0.1 0.2 0.4 0.8 ]},'txt', ...
                       {{num2str(roundsig(val(1),2)), '0.1','0.2','0.4','0.8'}});
else
    % score ticks are provided
    [AT1, BT1] = get_labels(opts.gsLabels);
    if opts.qv_thresh > 0
        score_ticks=struct('vals',{[val(1) AT1 ]},'txt',{{num2str(roundsig(val(1),2)), BT1{:}}});
    else
        score_ticks=struct('vals',{AT1},'txt',{BT1});
    end
end

%% SUBFUNCTION - generate tick structure for log log scale

function qv_ticks = loglog_q_ticks(minQ,nsubdiv,sigQ,spacing)

if ~exist('nsubdiv','var') || isempty(nsubdiv)
    nsubdiv = 6;
end
if ~exist('sigQ','var') || isempty(sigQ)
    sigQ = .25;
end
if ~exist('spacing','var') || isempty(spacing)
    spacing = .10;
end

%% find "nice" tick values for log log scale

% define scaling from q-value to plot
qscale = @(q) log10(-log10(q)+1);

span = qscale(minQ);
ruffpos = span .* ((1:nsubdiv)/nsubdiv);
expos = 10 .^ (ruffpos)-1;

% round exponents to one significant digit
logexpos_exp = floor(log10(expos));
logexpos_mant = 10 .^ (log10(expos)-logexpos_exp);
scale1 = [round(logexpos_mant(1:end-1)) ceil(logexpos_mant(end))];
scale1 = scale1 .* 10 .^ logexpos_exp;

% eliminate values < 10^-1
scale1(scale1<1) = [];

% eliminate values that encroach on minimum Q
if ~isempty(scale1)
    pos1 = log10(scale1+1);
    mindist = pos1(end) * spacing;
    scale1(abs(pos1-qscale(sigQ))<mindist) = [];
end
% eliminate duplicate values
scale1(scale1==circshift(scale1,[0 -1])) = [];

% minimum 10^-1 for empty scales
if isempty(scale1)
    scale1 = 1;
end
    
% create output structure
scale_text = strcat('10^{-',cellstr(num2str(scale1(:))),'}');
scale_text = regexprep(scale_text,' ','');
tick_vals = 10.^(-scale1);

% prepend significance threshold
tick_vals = [sigQ tick_vals];
scale_text = [cellstr(num2str(sigQ));scale_text]';

% create tick structure
qv_ticks=struct('vals',{tick_vals},'txt', {scale_text});

