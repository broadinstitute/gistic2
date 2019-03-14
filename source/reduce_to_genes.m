function G = reduce_to_genes(D,rg,rg_field,options)
% REDUCE_TO_GENES Convert SNP copy number data to gene level data
%
%   G = reduce_to_genes(D,RG,RG_FIELD,options)
%
%   Converts SNP-level copy number data in D to gene level copy data in G
%   using the reference genome in structure array RG. RG_FIELD is a string
%   specifying the field in RG to use as a gene identifier, default 'symb'.
%   OPTIONS is a struct whose fields are optional parameters:
%       OPTIONS.collapse_method can be 'all' (default), 'mean', 'median',
%               'min', 'max', or 'extreme'
%       OPTIONS.find_snps_type is the type argument passed to find_snps
%               0 for SNP markers strictly contained in the gene region
%               1 to allow flanking SNPS to give gene CN data
%               2 to use flanking SNPS only if none are strictly contained
%       (if the genes in RG have a 'snps' entry, then 
%

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


    if ~exist('options','var')
        options = struct;
    end
    if ~isfield(options,'collapse_method')
        options.collapse_method = 'mean';
    end
    if ~isfield(options,'find_snps_type')
        options.find_snps_type = 1;
    end
    if ~exist('rg_field','var') || isempty(rg_field)
        rg_field = 'symb';
    end

    % compress genes with same identifier in reference genome 
    ids_chr = strcat({rg.(rg_field)},'|chr',num2chromosome([rg.chrn]));
    [uqids,ui,oi] = lunique(ids_chr);
    
    % leave '|chrN' for genes that occur on multiple chromosomes
    ids = {rg(ui).(rg_field)};
    unids = unique(ids);
    [~,x1,x2] = match_string_sets_hash(ids,unids);
    mm = sparse(x2,x1,1);
    dup_chr = full(sum(mm,2)) > 1;
    unique_chr = ~(dup_chr'*mm);
    uqids(unique_chr) = ids(unique_chr);
    
    
    [~,m1,m2] = match_string_sets_hash(uqids(oi),uqids); %!
    gmatch = sparse(m1,m2,true); %!
%!  [gmatch,m1,m2]=match_string_sets_hash(ids,uqids);
    % m1 indexes ids
    % m2 corresponds ids, indexes uqids
    
    % create mapping from compresssed reference genome to snps
    usize = [length(uqids) 1];
    snps = cell(usize);    % cell array of snp positions
    rgindex = cell(usize); % cell array of reference genome indices
    chrn = zeros(usize);   % chromosome (NaN if ambiguous)
    g_start = Inf(usize);  % minimum gene start
    g_end = zeros(usize);  % maximum gene end
    snp_start = Inf(usize);% minimum flanking snp start
    snp_end = zeros(usize);% maximum flanking gene end
    
    % loop over unique gene IDs
    for ag = 1:length(uqids)
        if mod(ag,1000)==0
            verbose('reduce_to_genes: defining gene %d of %d',30,ag,length(uqids));
        end
        matches = find(gmatch(:,ag));
        rgindex{ag} = matches;
        if ~isempty(matches)
            % loop over all matches to unique gene ID
            for k = matches'
                g_start(ag) = min(g_start(ag),rg(k).start);
                g_end(ag) = max(g_end(ag),rg(k).end);
                if chrn(ag) == 0
                    chrn(ag) = rg(k).chrn;
                elseif chrn(ag) ~= rg(k).chrn;
                    assert(false);
                    %!!! this should no longer happen since we give genes
 %!                 % on different chromosomes distinct identifiers
 %!                 verbose('Conflicting chromosomes (%d and %d) for gene %s !',30, ...
 %!                             chrn(ag),rg(k).chrn,rg(k).(rg_field));
 %!                 chrn(ag) = min(chrn(ag),rg(k).chrn);
                end
                if g_end(ag) > g_start(ag) && ~isnan(chrn(ag))
                    % combine markers for the same gene
                    if isfield(rg,'snps')
                        snps{ag} = union(snps{ag},rg(k).snps);
                    else
                        snps{ag} = union(snps{ag},find_snps(D,chrn(ag),...
                                     g_start(ag),g_end(ag),options.find_snps_type));
                    end
                else
                    verbose('Bad gene data for %s !',30, ...
                            rg(k).(rg_field));
                end
            end
            % get positions of (possibly flanking) snps
            if ~isempty(snps{ag})
                snp_start(ag) = D.pos(min(snps{ag}));
                snp_end(ag) = D.pos(max(snps{ag}));
            end
        end
    end
    
    % create index for compressing/ordering genes/snps
    geneHasSnps = find(cellfun(@(s) ~isempty(s),snps));
    snps = snps(geneHasSnps);
    firstsnps = cellfun(@(s) s(1),snps);
    [firstsnps,order] = sort(firstsnps);
    geneHasSnps = geneHasSnps(order);
    snps = snps(order);

    
    % create output struct by removing matrix fields
    matrix_fields = {'dat','orig','cbs','affy_calls'}; % TODO add more?
    if isfield(D,'matrix_fields')
        matrix_fields = unique([matrix_fields D.matrix_fields]);
    end
    G = rmfield_if_exists(D,matrix_fields);
    % remove some other fields we won't need
    G = rmfield_if_exists(G,{'marker','chr','pos','cM', ...
                            'score','cbs_rl','orig', ...
                            'history','origidx','gorigidx',...
                            'medians','Qs'});
    % collapse G rows from snp to gene size 
    G = reorder_D_rows(G,firstsnps);
    % add some gene fields
    G = add_D_field(G,'gene',{'gdesc','snps','rgindex','chrn',...
                               'gstart','gend'});
    G.gdesc = uqids(geneHasSnps)'; %! achilles: geneID
    G.snps = snps;
    G.rgindex = rgindex(geneHasSnps);
    G.chrn = chrn(geneHasSnps);
    G.gstart = g_start(geneHasSnps); %! achilles: gene_start
    G.gend = g_end(geneHasSnps);     %! achilles: gene_end
    G.collapse_settings = options;

%{
    % create the big gene+location text description
    losnptxt = cellfun(@(x) num2str(min(x)),G.snps,'UniformOutput',false);
    hisnptxt = cellfun(@(x) num2str(max(x)),G.snps,'UniformOutput',false);
    G.gacc = strcat(G.gdesc,' - ',num2chromosome(G.chrn)',':',...
                    cellstr(num2str(G.gstart)),'-',cellstr(num2str(G.gend)),...
                    '[',losnptxt,'-',hisnptxt,']');
%}

    %% collapse copy number data

    % loop across all copy number matrix fields in D
    for field = matrix_fields
        fld = char(field);
        if isfield(D,fld)
            % add empty matrix field
            G.(fld) = nan(length(snps),size(D.(fld),2));
            % loop over genes
            for i = 1:length(snps)
                if mod(i,1000)==0
                    verbose('reduce_to_genes: reducing gene %d of %d',30,i,length(snps));
                end
                % compress markers to current gene
                switch options.collapse_method
                    case 'mean'
                        G.(fld)(i,:) = nanmean(D.(fld)(snps{i},:),1);
                    case 'median'
                        G.(fld)(i,:) = nanmedian(D.(fld)(snps{i},:),1);
                    case 'min'
                        G.(fld)(i,:) = nanmin(D.(fld)(snps{i},:),[],1);
                    case 'max'
                        G.(fld)(i,:) = nanmax(D.(fld)(snps{i},:),[],1);
                    case 'extreme'
                        min_val = nanmin(D.(fld)(snps{i},:),[],1);
                        max_val = nanmax(D.(fld)(snps{i},:),[],1);
                        vals = [min_val; max_val];
                        [~,mi] = nanmax(abs(vals));
                        v = zeros(1,length(min_val));
                        v(mi==1) = min_val(mi==1);
                        v(mi==2) = max_val(mi==2);
                        G.(fld)(i,:) = v;
                end
            end % loop over genes
        end
    end % loop over matrix fields

