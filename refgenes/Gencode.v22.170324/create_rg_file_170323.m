% 2017-03-23 create refgene for GENCODE v22 genes
% (derived from v18/create_test_rg_file_140127.m)


clear

% create test rg file
GG = read_R_table('gencode_genes.tsv');     % GENCODE genes
length(GG); %% 60483


%% gene filtering

GP = GG;

% no types IG_*, TR_*, antisense or sense_intronic
GP(strmatch('IG_',{GP.gene_type})) = [];
GP(strmatch('TR_',{GP.gene_type})) = [];
GP(ismember({GP.gene_type},{'antisense','sense_intronic'})) = [];

% no level 3 types except status=KNOWN miRNA, snoRNA, snRNA, scaRNA
lev3_except = strcmp({GP.gene_status},'KNOWN') & ismember({GP.gene_type},{'miRNA','snoRNA','snRNA','scaRNA'});
GP([GP.level]==3 & ~lev3_except) = [];

% remove lincRNAs unless status=KNOWN or level 1 status=NOVEL
GP(strcmp({GP.gene_type},'lincRNA') & ~(([GP.level] > 3 & strcmp({GP.gene_status},'KNOWN')) | ...
                                       ([GP.level] == 1 & ismember({GP.gene_status},{'KNOWN','NOVEL'})))) = [];

length(GP) % 40521
length(unique({GP.gene_name})) % 40380
                                   

%% deal with multiple genes with same gene_name
gene_names = {GP.gene_name};
unames = unique(gene_names);
[~,x1,x2] = match_string_sets_hash(gene_names,unames);
M = sparse(x1,x2,true);


%% deal with duplicate gene names
samenamex = find(full(sum(M,1)) > 1); % index of genes with duplicate names
ensembl_id = regexprep({GP.gene_id},'\.[0-9]*$','');
retain = true(size(ensembl_id));

% loop over names with multiple genes
for i = 1:length(samenamex)
    genes = find(M(:,samenamex(i))); % indices of matches to current duplicated name
    % single-number genomic coordinates
    starts = 1e10*chromosome2num({GP(genes).chr})+[GP(genes).start];
    ends = 1e10*chromosome2num({GP(genes).chr})+[GP(genes).end];
    % test each gene associated with the name for overlap
    for j = 1:length(genes)
        if retain(genes(j))
            overlaps = retain(genes) & (~(starts(j) > ends | ends(j) < starts));
            if sum(overlaps) > 1
                % overlapping genes: keep the largest, mark the rest for removal
                sizes = ends - starts;
                [~,mi] = max(sizes(overlaps));
                retain(genes(overlaps)) = false;
                retain(genes(mi)) = true;
            end
            % test (1) if gene is retained and (2) there is one other non-overlapping retained gene
            if retain(genes(j)) && sum(retain(genes) & ~overlaps) > 0
                 % non-overlapping, retained genes must be distinguished with unique suffix
                if sum(strcmp(GP(genes(j)).chr,{GP(genes).chr})) == 1
                    suffiq = GP(genes(j)).chr;
                else % (NOTE: not logically airtight, but works for current Gencode definitions)             
                    suffiq = ensembl_id{genes(j)}((end-5):end);
                end
                gene_names{genes(j)} = [GP(genes(j)).gene_name,'|',suffiq];
            end
        end
    end
end

% filter genes to remove overlapping duplicates
GP = GP(retain);
gene_names = gene_names(retain);
ensembl_id = ensembl_id(retain);

% NOTE: punting on the locus-link mapping because one locus_link may map to multiple ensembl IDs
%{
%% map locus_ids
N2E = read_R_table('ncbi2ensembl.tsv');     % NCBI-to-ENSMBL mappings
[~,nx,gx] = match_string_sets_hash({N2E.ensembl_id},ensembl_id);
locus_id = zeros(length(GP),1);
locus_id(gx) = [N2E(nx).locus_id];
locus_id(locus_id==0) = -cellfun(@(x) str2num(x(end-6:end)),ensembl_id(locus_id==0));
%}

% instead, we'll use last 7 digits of the ensembl ID
locus_id = cellfun(@(x) str2num(x(end-6:end)),ensembl_id);

%% package it into a struct
% NOTE: required fields are 'symb','chrn','start','end', and 'locus_id
% NOTE: some fields are cast to int32 to be compatible with GISTIC
rg = struct('symb',gene_names,...
            'chrn',num2cell(chromosome2num({GP.chr})),...
            'start',num2cell(int32([GP.start])),...
            'end',num2cell(int32([GP.end])),...
            'locus_id',num2cell(int32(locus_id)),...
            'gene_id',{GP.gene_id},...
            'gene_type',{GP.gene_type},...
            'gene_status',{GP.gene_status},...
            'source',{GP.source},...
            'chr',{GP.chr},...
            'strand',num2cell(double(strcmp({GP.strand},'+')))...
);

%% load cyto file
% use hg38 cytoband map from UCSC 
cyto = read_cytoband_data('cytoBand.txt'); % copy file locally
                                           %!cyto = read_cytoband_data('/xchip/gistic/variables//20160920_UCSC_dump/hg38/cytoBand.txt');

%!!! looks like hg38 gets overwritten by hg19 cyto file in 2017 code???
cyto_file = './hg19/cyto_ucsc_111024/cyto_hg19_111024.mat';
%!cyto_file ='/xchip/gistic/variables/cyto_ucsc_111024/hg19/cyto_hg19_111024.mat';

%!load(cyto_file); %!!! this code was near the end

%% create metadata

rg_path = '/xchip/beroukhimlab/gistic/reference_genomes/gencode/v22/';
%!rg_path = '/xchip/gistic/variables/gencode/v22/';
output_file = ['GENCODE.v22.',date6,'.refgene.mat'];
make_ver = '0.9 (beta)';

url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz';

[~,me] = unix('echo $USER');

% metadata structure
rg_info = struct( ...
    'assembly',    'GENCODE.v22', ...
    'source',      'GENCODE', ...
    'path',         rg_path, ...
    'url',          url,...
    'script',       [mfilename('fullpath') '.m'], ...
    'make_version', make_ver, ...
    'maker',        strtrim(me), ...
    'date',         date6 ...
);

% create output directory if necessary
if ~exist(rg_path)
    mkdir(rg_path)
end

%% save file
save([rg_path output_file],'-v7.3','rg','rg_info','cyto');
