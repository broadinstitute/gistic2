function [D,nsnpsremoved,m2] = remove_cnv(D,cnv_file,byposition)
%REMOVE_CNV remove germline copy number variants from data
%
%    [D,NREMOVED,REMIDX] = REMOVE_CNV(D,CNV_FILE,BYPOSITION)
%
% Remove markers in D_struct D that are contained by segments with high
% copy number variation specified in the file named by CNV_FILE. The
% returned D_struct has the markers covered by the specified CNVs removed
% alond with the associated data. The optional NREMOVED and REMIDX outputs
% indicate the number of markers removed and their index within the input
% markers.
%
% Two CNV_FILE formats will work with this function. The "by position"
% format has a 6-columns: (1) ID, (2) chromosome, (3) base position start,
% (4) base position end, (5) flanking start, and (6) flanking end. The
% information from columns 2, 5, and 6 is used to determine genomic
% segments to exclude.
%
% The non-positional CNV file format has marker identifiers in the first
% column that are matched against the markers in the D struct. This
% unsegmented format is specific to the hybridization platform.
%
% The optional BYPOSITION argument forces the file interpretation, 1= by
% position, 0 = by marker. If onitted, the format is determined by the
% number of columns in the file.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% iterator if cnv_file is a cell array
if iscell(cnv_file)
    if ~exist('byposition','var')
        byposition = [];
    end
    m2 = [];
    for i=1:numel(cnv_file)
        [D,~,idx] = remove_cnv(D,cnv_file{i},byposition);
        m2 = union(m2,idx);
    end
    nsnpsremoved = length(m2);
    return
end

% make sure the cnv file exists
if ~exist(cnv_file,'file')
    throw(MException('snp:remove_cnv:FileNotFound',...
                     'CNV file ''%s'' not found.',cnv_file));
end
% confirm we can open the CNV file
fid = fopen(cnv_file);
if fid < 0
    throw(MException('snp:remove_cnv:cnvfileOpenError',...
                     'Cannot open CNV file ''%s''.',cnv_file));
end
% read the first line
str = fgets(fid);
firstline = textscan(str,'%s%s%s%s%s%s','Delimiter',char(9));
fclose(fid);

% test for PC/unix text file (line delimited by '\r\n' or '\n' versus Mac '\r')
pcunix_text = str(end) == char(10);

if ~exist('byposition','var') || isempty(byposition)
    % detect positional format if at least 6 columns
    byposition = true;
    if isempty(firstline{6}) || isempty(firstline{6}{1})
        byposition = false;
    end
end

%% test for header
if byposition
    % assume header if third column of position-based file is not numeric
    hasheader = isnan(str2double(firstline{3}));
else
    % assume header if the first row/column is not among the markers
    hasheader = isfield(D,'marker') && ~any(strcmp(firstline{1},D.marker));
end

if byposition
    %% read positional CNV file data
    % data columns: (1) ID, (2) chromosome, (3) start position,
    % (4) end position, (5) flanking start, (6) flanking end
    fid = fopen(cnv_file);
    if hasheader
        fgetl(fid);
    end
    % scan in the CNV data columns
    % (NOTES: only the first 6 columns are read, remaining text to EOL
    %  is skipped.)
    try
        % choose format string based on line termination convention
        if pcunix_text
            fmt = '%s%s%n%n%n%n%*[^\n]';
        else
            fmt = '%s%s%n%n%n%n%*[^\r]';
        end
        cnv = textscan(fid,fmt,'Delimiter',char(9),...
                      'ReturnOnError',0,...'CommentStyle',char(13),...
                      'TreatAsEmpty','NA','EmptyValue',NaN);
    catch me
        fclose(fid);
        throw(MException('snp:remove_cnv:SegmentInputScanError',...
                    strcat('Error scanning segment file ''%s'':\n',me.message),...
                    cnv_file));
    end
    fclose(fid);

    % test for missing flanking start or end positions
    if any(isnan(cnv{5})) || any(isnan(cnv{6}))
        badline = find(isnan(cnv{5})|isnan(cnv{6}),1,'first') + hasheader;
        throw(MException('snp:remove_cnv:missingtPositionData',...
                'Data missing from line %d of CNV file ''%s''',badline,cnv_file));
    end

    cnvchrn = chromosome2num(cnv{2});
    % if errors, generate message about the first error
    if any(isnan(cnvchrn))
        badline = find(isnan(cnvchrn),1,'first');
        throw(MException('snp:remove_cnv:BadSegment',...
                ['%d invalid chromosome(s) detected in segment file ''%s''.\n',...
                'First invalid chromosome is ''%s'' on line %d'],...
                sum(isnan(cnvchrn)),cnv_file,cnv{2}{badline}, badline+hasheader));
    end
    
    % CNVs transformed into ordered genomic space using flanking start/end
    cnvstartpos = cnvchrn*1e11 + double(cnv{5});
    cnvendpos = cnvchrn*1e11 + double(cnv{6});
    % flip any reversed segments
    begends = [cnvstartpos,cnvendpos];
    [begends,colord] = sort(begends,2);
    flipped = colord(:,1)==2;
    if any(flipped)
        % warn about reversed segments
        firstflip = find(flipped,1,'first') + hasheader;
        warning( 'snp:remove_cnv:ReverseSegment',...
                ['%d segments are reversed in CNV file ''%s''.\n' ...
                 'The first reversed segment is on line %d.'], ...
                sum(flipped), cnv_file, firstflip);
    end
    cnvstartpos = begends(:,1);
    cnvendpos = begends(:,2);

    %% identify markers that fall within CNV segments
    
    % sort CNVs by start position
    [cnvstartpos,si] = sort(cnvstartpos);
    cnvendpos = cnvendpos(si);
    % calculate marker positions
    markerpos = double(D.chrn)*1e11 + double(D.pos);
    % define values used as reference counter increments
    % for entering or leaving a CNV segment
    cnvstartval = ones(size(cnvstartpos));
    cnvendval = -ones(size(cnvendpos));
    markerval = zeros(size(markerpos));
    % create aligned positions and reference increments
    allpos = [cnvstartpos; markerpos; cnvendpos];
    allval = [cnvstartval; markerval; cnvendval];
    % sort the reference increments according to position
    [~,posidx] = sort(allpos);
    allval = allval(posidx);
    % sum the increments to CNV reference count
    cnvcount = cumsum(allval);
    % marker reference counts > 0 are in CNVs
    cnvmarkers = cnvcount(allval==0) > 0;
    
    %% remove CNV markers

    % if C.dat is a SegArray, use logical SegArray to represent the CNVs
    use_segarray = isfield(D,'dat') && isa(D.dat,'SegArray');
    if use_segarray
        cnvmarkers = SegArray(cnvmarkers);
    end
    % remove the data
    verbose(['Removing ' num2str(sum(cnvmarkers)) ' markers'],10);
    D=reorder_D_rows(D,~cnvmarkers);
    % third output is numeric index
    m2 = find(cnvmarkers);

else
    %% non-positional CNV file
    if isfield(D,'marker')
        fid = fopen(cnv_file);
        if hasheader
            fgetl(fid);
        end
        % first column is of marker names
        cnv=textscan(fid,'%s%*[^\n\r]','ReturnOnError',0);
        fclose(fid);
        % match markers
        [~,m1,m2]=intersect(cnv{1},D.marker);
        markers_not_found = length(cnv{1}) - length(m2);
        if logical(markers_not_found)
            missed = setdiff(1:length(cnv{1}),m1);
            warning('snp:remove_cnv:MissedMarkers', ...
                    ['%d of the %d markers listed in the CNV file ''%s''' ...
                     ' were not found in the dataset.\n',...
                     'The first missing marker is on line %d.'],...
                    markers_not_found,length(cnv{1}),cnv_file,missed(1)+hasheader);
        end
        verbose('Removing %d markers',10,length(m2));
        D = reorder_D_rows(D,setdiff(1:size(D.dat,1),m2));
    else
        warning('snp:remove_cnv:no_data_markers', ...
                ['The CNV file ''%s'' is marker-based, however, the ' ...
                 'data do not use markers! CNV file ignored.'],cnv_file);
        m2 = [];
    end
end

%% 2nd output derived from 3rd
nsnpsremoved = length(m2);

