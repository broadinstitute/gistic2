function refgene = load_refgene(refgene)
%LOAD_REFGENE - load reference genome and associated info

% GISTIC software version 2.0
% Copyright (c) 2011,2014 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

global REFGENE_INFO

if ischar(refgene)
    if exist(refgene,'file')
        load(refgene);
    else
        throw(MException('snp:load_refgene:BadFile','refgene file ''%s'' not found.'));
    end

    % ensure refgene file is valid
    if ~exist('rg','var') || ~exist('cyto','var')
        throw(MException('snp:load_refgene:BadFile','''%s'' is not a valid refgene file.'));
    end

    % if no RefGeneInfo, object, make up an hg refgene 
    if ~exist('RGI','var') || ~isa(RGI,'RefGeneInfo')
        RGI = RefGeneInfo;
        % newer refgenes have information about themselves
        if exist('rg_info','var')
            RGI = add_info(RGI,rg_info);
        end
    %!else
    %!    REFGENE_INFO = RGI;
    end
    
    % patch: some old refgenes assignes unmapped contigsnumbers
    bad_ref_genes = ~([rg.chrn] <= RGI.nchr); % (nans too)
    if any(bad_ref_genes)
        verbose('Removing bad reference genes.',30);
        rg(bad_ref_genes) = [];
    end

    % create structure for output
    refgene = struct;
    refgene.gene = rg;
    refgene.cyto = cyto;
    refgene.rgi = RGI;
    
    REFGENE_INFO = RGI;
elseif isstruct(refgene) && all(isfield(refgene,{'gene','cyto','rgi'}));
    % already loaded structure, just set global refgene info
    REFGENE_INFO = refgene.rgi;
else
    throw(MException('snp:load_refgene:Refgene','Invalid reference genome.'));
end

