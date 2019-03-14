function tab = genetables(D,A,~,rg,cyto,regs,calls,~,ext,...
                          ~,base_dir,genepattern,partial_hits,fname)
%GENETABLES output amp_genes, del_genes, table_amp and table_del tables

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~exist('base_dir','var')
    base_dir = [];
end

if ~exist('ext','var') 
    ext = '';
end
if ~exist('fname','var')
    fname = '';
end
  
if mean(D.pos)<10000 % is in Mb
    error('DOES NOT WORK WITH MB ANYMORE');
end

if ~isempty(A) && ischar(A.gacc)
    A.gacc=cellstr(A.gacc);
end

if ~exist('partial_hits','var') || isempty(partial_hits)
  partial_hits = ones(length(calls),1);
else
  if length(partial_hits) ~= length(calls)
    if length(partial_hits)== 1
      partial_hits = repmat(partial_hits,length(calls),1);
    else
      error(['Length of partial_hits does not match number of ' ...
             'amplifications/deletions/LOH categories tested'])
    end
  end
end

ampdel = {'amp','del'};

%% create data to output to table_{amp|del}.txt output file

tab=cell(1,length(calls));
for k=1:length(calls)  %calls is cell with 1/0 amplifications{1} and deletions{2} data
  tab{k}=[];
  for pki=1:size(calls{k},1)  % :-1:1  %loop over peaks
    tab{k}.idx(pki)=pki;
    %regs has all the info about the peaks
    %Add enl_peak_st field to regs ("enlarged peak region"?)
    if isfield(regs{k}(pki),'peak_wide_st')
      regs{k}(pki).enl_peak_st=regs{k}(pki).peak_wide_st;  
      regs{k}(pki).enl_peak_en=regs{k}(pki).peak_wide_en;
    else
      [regs{k}(pki).enl_peak_st,regs{k}(pki).enl_peak_en]=...
          enlarge_range(D,regs{k}(pki).peak_st,regs{k}(pki).peak_en,15,15,...
                        regs{k}(pki).st,regs{k}(pki).en);
    end
    
    [~,~,enl_region_to_next_snp] = genomic_location(D,{[regs{k}(pki).enl_peak_st ...
                        regs{k}(pki).enl_peak_en]},cyto,1);
    [~,~,peak_region_to_next_snp] = genomic_location(D,{[regs{k}(pki).peak_st ...
                        regs{k}(pki).peak_en]},cyto,1);
    % get unique genes in region
    [genes_in_region,closest_to_region] = genes_at(rg,D.chrn(regs{k}(pki).enl_peak_st), ...
                             enl_region_to_next_snp{1}(1), enl_region_to_next_snp{1}(2),1,partial_hits(k));
    genes_in_region = filt_uniq_gene(rg, genes_in_region);
    
    % get unique genes in peak
    [genes_in_peak,closest_to_peak] = genes_at(rg,D.chrn(regs{k}(pki).peak_st), ...
                           peak_region_to_next_snp{1}(1),peak_region_to_next_snp{1}(2),1,partial_hits(k));
    genes_in_peak = filt_uniq_gene(rg, genes_in_peak);
        
    tab{k}.chr(pki)=D.chrn(regs{k}(pki).peak_st);
    tab{k}.region_start(pki)=D.pos(regs{k}(pki).st);
    tab{k}.region_end(pki)=D.pos(regs{k}(pki).en);
    
    tab{k}.peak_start(pki)=D.pos(regs{k}(pki).peak_st);
    tab{k}.peak_end(pki)=D.pos(regs{k}(pki).peak_en);
    
    tab{k}.enlarged_peak_start(pki)=D.pos(regs{k}(pki).enl_peak_st);
    tab{k}.enlarged_peak_end(pki)=D.pos(regs{k}(pki).enl_peak_en);
    
    tab{k}.n_genes_in_region(pki)=length(genes_in_region);
    tab{k}.n_genes_in_peak(pki)=length(genes_in_peak);
    
    if ~isempty(genes_in_region)
      tab{k}.genes_in_region{pki}=sprintf('%s,',rg(genes_in_region).symb);
    elseif ~isempty(closest_to_region)
      tab{k}.genes_in_region{pki}=['[' rg(closest_to_region).symbol ']'];
    else
      tab{k}.genes_in_region{pki}='Empty';
    end
    
    if ~isempty(genes_in_peak)
      tab{k}.genes_in_peak{pki}=sprintf('%s,',rg(genes_in_peak).symb);
    elseif ~isempty(closest_to_peak)
      tab{k}.genes_in_peak{pki}=['[' rg(closest_to_peak).symbol ']'];
    else
      tab{k}.genes_in_peak{pki}='Empty';
    end    
  end
end


%% output the table_{amp|del}.txt file

if ~exist('genepattern','var') || isempty(genepattern) || ~genepattern
  for k=1:length(calls)
    f = fopen([base_dir fname 'table_' ampdel{k}  ext '.txt'],'w');
    fprintf(f,['index\tchromosome\tregion_start\tregion_end\tpeak_start\tpeak_end\tenlarged_peak_start\tenlarged_peak_end\t' ...
        'n_genes_in_region\tgenes_in_region\tn_genes_in_peak\tgenes_in_peak\tn_genes_on_chip\' ...
        'tgenes_on_chip\ttop 3\n']);
    if ~isempty(tab{k}) && isfield(tab{k},'idx')
        for i=1:length(tab{k}.idx)
            if ~isfield(tab{k},'chr')
                fprintf(f,'%d',tab{k}.idx);
            else
                fprintf(f,['%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t' ...
                    '%d\t%s\t%d\t%s\t%d\t%s\t%s'],...
                    tab{k}.idx(i),...
                    num2chromosome(tab{k}.chr(i)),...
                    tab{k}.region_start(i),...
                    tab{k}.region_end(i),...
                    tab{k}.peak_start(i),...
                    tab{k}.peak_end(i),...
                    tab{k}.enlarged_peak_start(i),...
                    tab{k}.enlarged_peak_end(i),...
                    tab{k}.n_genes_in_region(i),...
                    tab{k}.genes_in_region{i},...
                    tab{k}.n_genes_in_peak(i),...
                    tab{k}.genes_in_peak{i});
%                     tab{k}.n_genes_on_chip(i),...
%                     tab{k}.genes_on_chip{i},...
%                     tab{k}.top3{i});
            end
            if isfield(tab{k},'md')
                for jj=1:size(tab{k}.md,2)
                    fprintf(f,'\t%s',tab{k}.md{i,jj});
                end
            end
            fprintf(f,'\n');
        end
    end
    fclose(f);
  end
end

%%  write {amp|del}_genes.txt outputs

% create sort order by residual q-value
writeidx = {[],[]};
for k = 1:length(regs)
  if ~isempty(regs{k})
    [~,writeidx{k}] = sort([regs{k}.resid_qv]);
%! no longer by genomic location
%!  [~,writeidx{k}] = sort(1e11.*double([regs{k}.chrn])+double([regs{k}.peak_st]));
  end
end



%change row or label here.
%must have an entry in rows and labels for eachfield. 
%if fields are added, need to also modify 'make alltext cell' section below
rows.cy = 1;
labels.cy = 'cytoband';
rows.qv = 2;
labels.qv = 'q value';
rows.rqv = 3;
labels.rqv = 'residual q value';
rows.wp = 4;
labels.wp = 'wide peak boundaries';
rows.gns = 5;
labels.gns = 'genes in wide peak';

for k = 1:length(regs)  %Do for amps and dels
  % open gene file
  filename = [base_dir fname ampdel{k} '_genes' ext '.txt'];
  fid = fopen(filename,'w');
    
  if ~isempty(regs{k})
    genesinregion = cellfun(@(x) regexp(x,'[^,]+','match'),tab{k}.genes_in_region,'UniformOutput',0);
    
    % get cytoband names
    lesionlocs = double(D.pos([regs{k}.peak])') + 1e11.*double([regs{k}.chrn]);  %where the lesions start
    cytolocs = double([cyto.start]) + double([cyto.chrn].*1e11);  %where the cytobands start
    [srtcytolocs,ci] = sort(cytolocs);
    [srtll,li] = sort(lesionlocs);  %indices needed to get lesions in order
    [~,si] = sort([srtcytolocs srtll]);  %
    cytobands = {cyto(ci(find(si>length(cytolocs))-(1:length(lesionlocs)))).name};
    cytobands(li) = cytobands;
  
    % make alltext output cell array
    alltext = repmat({''},length(fieldnames(rows))+max(cellfun(@length,genesinregion))-1,length(cytobands)+1);
    alltext(cell2mat(struct2cell(orderfields(rows))),1) = struct2cell(orderfields(labels));
    alltext(rows.cy,2:end) = cytobands;
    alltext(rows.qv,2:end) = cellfun(@num2str,{regs{k}.qv},'UniformOutput',0);
    alltext(rows.rqv,2:end) = cellfun(@num2str,{regs{k}.resid_qv},'UniformOutput',0);
    tmp_cell=genomic_location(D,mat2cell([cat(1,regs{k}.peak_wide_st) ...
                        cat(1,regs{k}.peak_wide_en)],ones(length(regs{k}),1),2),cyto,1);
    if ~iscell(tmp_cell)
      tmp_cell={tmp_cell};
    end
    alltext(rows.wp,2:end) =  tmp_cell;
    alltext(rows.wp,2:end) = regexprep(alltext(rows.wp,2:end),'\(.+\)','');
    
    
    filledgenesinregion = vertcatfill(genesinregion{:});
    filledgenesinregion = reshape([filledgenesinregion{:}],length(filledgenesinregion{1}),length(filledgenesinregion));
    alltext(rows.gns:end,2:end) = filledgenesinregion;

    % reorder output and write file
    alltext = alltext(:,[1 writeidx{k}+1]);
    alltextvect = reshape(alltext',1,numel(alltext));
    fprintf(fid,[repmat('%s\t',1,size(alltext,2)) '\n'],alltextvect{:});
  else
    % no peaks - emit empty gene file
    fprintf(fid,'%s\n%s\n%s\n%s\n%s\n',labels.cy,labels.qv,labels.rqv,labels.wp,labels.gns);
  end
  % close gene file
  fclose(fid);
end

%% local functions

% filter peak gene selection to include unique genes
% Input: rg is reference genome, gene_select is a vector of indexing the
% rg genes in the peak
% Returns: filtered indices, optionally indices used to filter the indices
function [gene_select,pique] = filt_uniq_gene(rg,gene_select)
if ~isempty(gene_select)
    if isfield(rg,'gene_id')
      % if available, use gene_id string identifier to partition
      [~,pique,~] = lunique({rg(gene_select).gene_id}');
    elseif isfield(rg,'locus_id')
      % otherwise use legacy locus_id numeric identifier
      [~,pique,~] = lunique([rg(gene_select).locus_id]');
    else
      % fallback is to use symbol
      [~,pique,~] = lunique({rg(gene_select).symb}');
    end
    gene_select = gene_select(pique);
end
