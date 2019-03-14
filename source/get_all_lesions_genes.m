function gene_names = get_all_lesions_genes(all_lesions_file,regs,reg_names)
%GET_ALL_LESIONS_GENES get gene_names from the all_lesions file

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~isempty(all_lesions_file) && exist(all_lesions_file,'file') 
  %gene names are supplied, so the gene names are used
    fid = fopen(all_lesions_file);
    line = fgetl(fid);
    header = textscan(line,'%s','Delimiter','\t');
    
    % detect the version of the all_lesions file
    if strcmp(header{1}{1},'Lesions')
        oldversion = 1;
        fclose(fid);
    elseif strcmp(header{1}{1},'Unique Name')
        headerlength = length(header{:});
        oldversion = 0;
        M = textscan(fid,repmat('%s',1,headerlength),'Delimiter','\t');
        fclose(fid);
    else
        error('Unrecognized format for all lesions file')
    end
    
    
    if ~oldversion
        gene_names = {[] []};
        descriptoridx = strmatch('Descriptor',header{1},'exact');
   
        allidxs = find(~cellfun(@isempty,regexp(M{strmatch('Unique Name',header{1},'exact')},'[\d]$')));
        ampidx = allidxs(strmatch('Amplification Peak',M{strmatch('Unique Name',header{1},'exact')}(allidxs)));
        delidx = allidxs(strmatch('Deletion Peak',M{strmatch('Unique Name',header{1},'exact')}(allidxs)));
                      
        gene_names{1} = M{descriptoridx}(ampidx)';
        gene_names{2} = M{descriptoridx}(delidx)';
    else
        M = read_dlm_file(all_lesions_file, '\t'); %read in all lesions file

        if ~isempty(regs) && size(regs{1},2)>0
            %create a structure indicating where the amplified regions are in the all lesions file
            for k=1:size(M,2)
                interm_Amp = strmatch('Any-AP',M{k}{1});
                if interm_Amp
                    amp_Gene_IDX(k) = interm_Amp;
                end
            end

            for repeat = 1:2 %add repeat to initialize the geneListA structure the first time around
                for p=1:size(amp_Gene_IDX,2)
                    if amp_Gene_IDX(p) %if a region of amplification
                        for r = 1:size(reg_names{1},2)

                            %see where the cytoband is present in the all lesions file
                            found = findstr(reg_names{1}{r},M{p}{1});

                            if found >0 %if cytoband is present (i.e. has an index >0)
                                amp1 = M{p}{1};
                                amp2 = textscan(amp1,'%s','delimiter','\t');
                                amp3 = amp2{1}{2};
                                amp4 = textscan(amp3,'%s','delimiter',' ');
                                matchGenesA{1}{p} = amp4{1}{1};

                                %create structure with gene names in the location of amp_Gene_IDX
                                matchGenesA{1}{p} = matchGenesA{1}{p}(findstr(amp4{1}{1},'Any-')+4:end);

                                %create a field in the regs file with the gene name
                                regs{1}(r).genes = matchGenesA{1}{p};

                                %refine matchGenesA to get rid of empty locations
                                if ~isempty(matchGenesA{1}{p})
                                    geneListA{p-(size(matchGenesA{1},2)-size(regs{1},2))} = matchGenesA{1}{p};

                                    %refine geneListA to get rid of the | that is included in the all lesions file
                                    if geneListA{p-(size(matchGenesA{1},2)-size(regs{1},2))}(end) == '|'
                                        geneListA{p-(size(matchGenesA{1},2)-size(regs{1},2))} = ...
                                            geneListA{p-(size(matchGenesA{1},2)-size(regs{1},2))}(1:end-1);
                                    end

                                end
                            end
                        end
                    end
                end
            end

            geneListA = geneListA(amp_IDX); %sort out the gene list based on how regs was sorted
            gene_names{1} = geneListA;
        else
            gene_names{1} = [];
        end

        if ~isempty(regs) && size(regs{2},2)>0
            %create a structure indicating where the deleted regions are in the all lesions file
            for k=1:size(M,2)
                interm_Del = strmatch('Any-DP',M{k}{1});
                if interm_Del
                    del_Gene_IDX(k) = interm_Del;
                end
            end

            for repeat = 1:2 %add repeat to initialize the geneListA structure the first time around
                for p=1:size(del_Gene_IDX,2)
                    if del_Gene_IDX(p)
                        for r = 1:size(reg_names{2},2)

                            %see where the cytoband is present in the all lesions file
                            found = findstr(reg_names{2}{r},M{p}{1});
                            if found >0 %if cytoband is present (i.e. has an index >0)
                                del1 = M{p}{1};
                                del2 = textscan(del1,'%s','delimiter','\t');
                                del3 = del2{1}{2};
                                del4 = textscan(del3,'%s','delimiter',' ');
                                matchGenesD{1}{p} = del4{1}{1};

                                %create structure with gene names in the location of amp_Gene_IDX
                                matchGenesD{1}{p} = matchGenesD{1}{p}(findstr(del4{1}{1},'Any-')+4:end);

                                %create a field in the regs file with the gene name
                                regs{2}(r).genes = matchGenesD{1}{p};

                                %refine matchGenesA to get rid of empty locations
                                if ~isempty(matchGenesD{1}{p})
                                    geneListD{p-(size(matchGenesD{1},2)-size(regs{2},2))} ...
                                        = matchGenesD{1}{p};

                                    %refine geneListA to get rid of the | that is included in the all lesions file
                                    if geneListD{p-(size(matchGenesD{1},2)-size(regs{2},2))}(end) == '|'
                                        geneListD{p-(size(matchGenesD{1},2)-size(regs{2},2))} = ...
                                            geneListD{p-(size(matchGenesD{1},2)-size(regs{2},2))}(1:end-1);
                                    end

                                end
                            end
                        end
                    end
                end
            end

            geneListD = geneListD(del_IDX); %sort out the gene list based on how regs was sorted
            gene_names{2} = geneListD;
        else
            gene_names{2} = [];
        end

        if size(gene_names{1},2)==size(reg_names{1},2)
            for p=1:size(gene_names{1},2)
                if gene_names{1}{p}(1:2) == 'AP'
                    gene_names{1}{p} = reg_names{1}{p};
                end
            end
        else
            error('All lesion file does not contain the same number of regions as the regs file');
        end

        if size(gene_names{2},2)==size(reg_names{2},2)
            for p=1:size(gene_names{2},2)
                if gene_names{2}{p}(1:2) == 'DP'
                    gene_names{2}{p} = reg_names{2}{p};
                end
            end
        else
            error('All lesion file does not contain the same number of regions as the regs file');
        end
    end
else  %%% gene names are not supplied
    gene_names=reg_names;
end % if we have an all_lesions file
