function [segs,hasheader] = read_segfile(segfile)
%% ensure input files exist

if ~exist(segfile,'file')
    throw(MException('snp:read_segfile:FileNotFound',...
                     'Segment file ''%s'' not found.',segfile));
end

%% test segfile
fid = fopen(segfile);
if fid < 0
    throw(MException('snp:read_segfile:segfileOpenError',...
                     'Cannot open segmented data file ''%s''.',segfile));
end
% read the first line
str = fgets(fid);
fclose(fid);

% test for PC/unix text file (line delimited by '\r\n' or '\n' versus Mac '\r')
pcunix_text = str(end) == char(10);

%% Read segmented data

verbose('Reading Seg File ''%s''',10,segfile);

hasheader=check_if_has_header(segfile,3);
first_line=read_dlm_file(segfile,char(9),1);

% open segment file
fid = fopen(segfile);
if fid < 0
    throw(MException('snp:read_segfile:segfileOpenError',...
                'Error opening segmented data file ''%s''.',segfile));
end

% scan formatted text
try
    if pcunix_text % PC or Unix newline delimited lines
        format=['%q%q%n%n' repmat('%*q',1,length(first_line{1})-5) '%f\n'];
        segcol = textscan(fid,format,'ReturnOnError',0,'TreatAsEmpty','NA',...
                    'Headerlines',0+hasheader,'MultipleDelimsAsOne',1,...
                    'Delimiter',char(9),'WhiteSpace','\r');
    else % Macintosh carriage return delimited lines
        format=['%q%q%n%n' repmat('%*q',1,length(first_line{1})-5) '%f\r'];
        segcol = textscan(fid,format,'ReturnOnError',0,'TreatAsEmpty','NA',...
                    'Headerlines',0+hasheader,'MultipleDelimsAsOne',1,...
                    'Delimiter',char(9));
    end
catch me
    fclose(fid);
    throw(MException('snp:read_segfile:SegmentInputScanError',...
                    strcat('Error scanning segment file ''%s'':\n',me.message),...
                    segfile));
end
fclose(fid);

% put columns into readable variables
segs = struct;
segs.sample = deblank(segcol{1});
segs.chr = segcol{2};
segs.start = segcol{3};
segs.end = segcol{4};
segs.copyN = segcol{5};

% test for missing start or end positions
if any(isnan(segs.start)) || any(isnan(segs.end))
    badline = find(isnan(segs.start)|isnan(segs.end),1,'first') + hasheader;
    throw(MException('snp:read_segfile:missingtSegmentData',...
                    'Data missing in segment file ''%s'', line %d',segfile,badline));
end

% convert chromosome text to numeric representation
segs.chrn = chromosome2num(segs.chr);
% if errors, generate message about the first error
if any(isnan(segs.chrn))
    badline = find(isnan(segs.chrn),1,'first');
    throw(MException('snp:read_segfile:BadSegment',...
        ['%d invalid chromosomes detected in segment file ''%s''.\n',...
        'First invalid chromosome is ''%s'' on line %d'],...
        sum(isnan(segs.chrn)),segfile,segs.chr{badline},badline+hasheader));
end
