classdef RefGeneInfo
%RefGeneInfo - class for reference genome metadata
   
% GISTIC software version 2.0
% Copyright (c) 2011,2014 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    properties (SetAccess = 'public', GetAccess = 'public')
        txt2num = containers.Map;  % map names to indicbuildes
        chr = struct;
        nchr = 24;
        info = struct;
    end
    methods (Access = 'public',Static = true)
        %% constructors

        function RGI = RefGeneInfo()
            global REFGENE_INFO
            if isa(REFGENE_INFO,'RefGeneInfo')
                % return global refgene info
                RGI = REFGENE_INFO;
            else                
                % no global refgene info yet - create default human
                % define number to text arrays
                RGI.chr.symb = {'1','2','3','4','5','6','7','8','9','10',...
                                '11','12','13','14','15','16','17','18',...
                                '19','20','21','22','X','Y'};
                RGI.nchr = length(RGI.chr.symb);
                RGI.chr.text = strcat('chr',RGI.chr.symb);
                % load the text-to-number map
                RGI.txt2num = containers.Map;
                for i=1:length(RGI.chr.symb)
                    RGI.txt2num(RGI.chr.symb{i}) = i;
                    RGI.txt2num(RGI.chr.text{i}) = i;
                end
                RGI.txt2num('23') = 23;
                RGI.txt2num('24') = 24;
                RGI.chr.autosomal = (1:24)<23;
            end
        end
        
        % from build_info file (like mutsig)
        function RGI = fromBuildInfoFile(fname)
            global REFGENE_INFO
            RGI = RefGeneInfo;
            %!!!TODO check if over-writing existing REFGENE_INFO
            BI = read_R_table(fname);
            BI = BI(strcmp({BI.use},'D'));
            BI([BI.num]) = BI;
            RGI.nchr = length(BI);
            RGI.chr.symb = {BI.chromosome};
            RGI.chr.text = {BI.chromosome};
            RGI.txt2num = containers.Map;
            for c = 1:length(BI)
                RGI.txt2num(BI(c).chromosome) = c;
                aliases = regexp(BI(c).alias,',','split');
                for i = 1:length(aliases)
                    RGI.txt2num(aliases{i}) = c;
                    if i == 1
                        %!!! TODO resolve build_info.txt file convention
                        RGI.chr.symb{c} = aliases{1};
%!                      RGI.chr.text{c} = aliases{1};
                    end
                end
            end
            RGI.chr.autosomal = [BI.male]==2 & [BI.female]==2;
            REFGENE_INFO = RGI;
        end        
    end
    
    methods (Access = 'public',Static = false)
    
        
        %% informational methods
        % maximum chromosome number
        function num = maxChromosomeNum(RGI)
            num = length(RGI.chr.symb);
        end
        
        %% methods for chromosome information
        
        % this function shall have the same behavior as chromosome2num 
        function num = getChromosomeNum(RGI,text)
            if ischar(text)
                if isKey(RGI.txt2num,text)
                    num = RGI.txt2num(text);
                else
                    num = nan;
                end
            elseif iscellstr(text)
                % map chromosome strings to numbers with containers.Map 
                num = nan(size(text));
                
                keyed = find(isKey(RGI.txt2num,text));                
                valc = values(RGI.txt2num,text(keyed));
                for i=1:length(keyed)
                    num(keyed(i)) = valc{i};
                end
                %! loop is 10x faster than parallel approaches below when
                %! compiled
%!              num(keyed) = [valc{:}];
%!              num(keyed) = cell2mat(valc);
            else
                % error
                throw(MException('RefGeneInfo:getChromosomeNum:InvalidArgument',...
                    'RefGeneInfo::getChromosomeNum cannot process arguments of type ''%s''',...
                    class(text)));
            end
        end
        
        % this function shall have the same behavior as num2chromosome
        % (but handle errors better)
        function text = getChromosomeText(RGI,num)
            if all(num <= maxChromosomeNum(RGI))
                if isscalar(num)
                    text = RGI.chr.symb{num};
                else
                    text = RGI.chr.symb(num);
                end
            else
                % error
                throw(MException('RefGeneInfo:getChromosomeText:InvalidArgument',...
                    'Chromosome index excedes number of chromosomes in reference genome'));
            end
        end
        
        % flag autosomes
        function isauto = isAutosomal(RGI,chrn)
            isauto = false(size(chrn));
            inrange = (chrn > 0) & (chrn <= RGI.nchr); 
            isauto(inrange) = RGI.chr.autosomal(chrn(inrange));
        end
        
        % add non-chromosomal information about refgene
        function RGI = add_info(RGI,rg_info)
            RGI.info = rg_info;
        end
        
        
    end
end
