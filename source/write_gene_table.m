function write_gene_table(fname,G,cyto,rg)
%WRITE_GENE_TABLE save gene-level copy number data to file
%
%   WRITE_GENE_TABLE(FNAME,G,CYTO,RG)
%
% Output G (copy number reduced to gene level) to the file FNAME, 
% using CYTO and RG to add locus ID and cytoband annotation columns.
%

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

nsamples = size(G.dat,2);
ngenes = size(G.dat,1);

fid = fopen(fname,'w');

% print header
fprintf(fid,'%s\t%s\t%s','Gene Symbol','Gene ID','Cytoband');
for j=1:nsamples
    fprintf(fid,'\t%s',G.sdesc{j});
end
fprintf(fid,'\n');

% define format string for a row
fmt = ['%s\t%s\t%s',repmat('\t%1.3f',1,nsamples),'\n'];

% loop over genes (genomic dimension of G)
for j=1:ngenes
    modi(j,1000)
    band = find([cyto.chrn] == G.chrn(j) & [cyto.start] <= G.gstart(j),1,'last');
    % determine gene identifier
    if isfield(rg,'gene_id')
        % first try preferred string identifier
        gene_id = rg(G.rgindex{j}(1)).gene_id; % (if multiple, use ast)
    elseif isfield(rg,'locus_id')
        % next try legacy numeric locus_id
        gene_id = num2str(rg(G.rgindex{j}(1)).locus_id);
    else
        % fallback is nothing
        gene_id = '';
    end
    fprintf(fid,fmt,G.gdesc{j},gene_id,cyto(band).name,G.dat(j,:)); 
end  
fclose(fid);
