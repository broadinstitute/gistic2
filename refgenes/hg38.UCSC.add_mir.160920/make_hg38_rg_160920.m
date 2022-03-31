% build hg38 reference genome from 2016-09-20 snapshot of UCSC files

clear
dbstop if error

% input path for this script / UCSC download directory
UCSC_dump_path = [pwd,filesep];
% date that the UCSC files were downloaded
snapshot_date = '20-Sep-2016';
make_ver = '0.98';
% metadata
url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/';

[~,me] = unix('echo $USER');
rg_info = struct( ...
    'assembly',    'hg38', ...
    'source',      'UCSC', ...
    'path',         UCSC_dump_path, ...
    'url',          url,...
    'script',       [mfilename('fullpath') '.m'], ...
    'make_version', make_ver, ...
    'maker',        strtrim(me), ...
    'date',         snapshot_date ...
);

% output path
refseq_path = '/xchip/gistic/variables/hg38/';
output_basename = 'hg38.UCSC.add_miR.160920.refgene';

refgene_fname = [UCSC_dump_path 'refGene.txt'];
reflink_fname =[UCSC_dump_path 'refLink.txt'];
refseqstatus_fname = [UCSC_dump_path 'refSeqStatus.txt'];
wg_fname = [UCSC_dump_path 'wgRna.txt'];

hg_output_fname = [refseq_path 'hu_refGene_' output_basename '.txt'];
rg = make_hg_file(hg_output_fname,refgene_fname,reflink_fname, ...
                  refseqstatus_fname,wg_fname,1);
rg=add_chrn(rg);
% eliminate haplotypes
rg(isnan([rg.chrn])|[rg.chrn]>24) = [];
cyto=read_cytoband_data([UCSC_dump_path 'cytoBand.txt']);
save([refseq_path output_basename '.mat'],'rg','cyto','rg_info');
savestruct(rg,[refseq_path output_basename '.txt']);
