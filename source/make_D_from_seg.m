function D = make_D_from_seg(segfile,markersfile,options,use_segarray)
%MAKE_D_FROM_SEG make "D" data structure from segmented CN data and markers
%
%   D = make_D_from_seg(SEGFILE,MARKERSFILE,OPTIONS) 
%
%   D = make_D_from_seg(SEGFILE,MARKERSFILE,USE_OLD,USE_SEGARRAY)
%       (old form of arguments for backwards compatibility)
%
% SEGFILE is a path to a file (string) containing segmented data
%    in whitespace delimited columns in this order: 
%        sample, chromosome, start base,end base, num markers, log CN
%
% MARKERSFILE is either a path to a file containing marker
%    locations or a cell array of already loaded marker positions.   
%
% OPTIONS is a structure containing optional arguments:
%   OPTIONS.use_segarray is set to 1 to reduce the memory footprint 
%       of D.dat by using a SegArray object. The default is 0, which uses
%       a full array.
%   OPTIONS.compress_markers should be set to 1 to remove redundant markers 
%       and strictly order them by position before matching segmented data.
%       The default is 0.
%   OPTIONS.suppress_history creates an output D that does not keep a history
%       of operations or copy of previous values. The default is 1.
%   OPTIONS.chr_field adds D.chr, a cell array of chromosome strings for 
%       each SNP to D equivalent to D.chrn. The default is 1. Set to 0 to 
%       save memory.
%   OPTIONS.marker_field adds D.marker, a cell array of marker identifiers,
%       to D. The markers come from the first column of MARKERSFILE. The
%       default is 1. Set to 0 to save memory.
%   OPTIONS.islog indicates whether the data from the markers is log ratio
%       or absolute copy number. OPTIONS.islog is copied to D.islog. The
%       default is to auto-detect from the data.
%   OPTIONS.isMB indicates whether the marker positions are in bases 
%       (default 0) or megabase units (1).
%
% USE_OLD and USE_SEGARRAY are logical arguments for the old form of this
%   function, which is maintained for backward compatibility but deprecated.
%   USE_OLD will generate an error if has a nonzero value. USE_SEGARRAY is
%   the equivalent of OPTIONS.use_segarray, which supersedes it. 

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
%% process input parameters

% options / use_old flag
if exist('options','var')
    if ~isstruct(options)
        % third argument not a struct => obsolete use_old
        if options
            % never use_old
            error('GISTIC:Unsupported',...
                ['The third ''use_old'' parameter for make_D_from_seg',...
                ' is no longer supported.']);
        else
            options = struct;
        end
    end
else
    options = struct;
end

% put old use_segarray argument into struct
if exist('use_segarray','var') && ~isfield(options,'use_segarray')
    options.use_segarray = use_segarray;
end
% provide defaults for undefined options
options = impose_default_value(options,'use_segarray',true);
options = impose_default_value(options,'compress_markers',false);
options = impose_default_value(options,'suppress_history',true);
options = impose_default_value(options,'chr_field',true);
options = impose_default_value(options,'marker_field',true);
%!options = impose_default_value(options,'islog',true); %!!!detect from data
options = impose_default_value(options,'isMB',false);

%% ensure input files exist
if ~exist(segfile,'file')
    throw(MException('snp:make_D_from_seg:FileNotFound',...
                     'Segment file ''%s'' not found.',segfile));
end
if ischar(markersfile) && ~exist(markersfile,'file')
    throw(MException('snp:make_D_from_seg:FileNotFound',...
                     'Marker file ''%s'' not found.',markersfile));
end

%% read markers file
if ~iscell(markersfile)
    % markersfile is a file name
    verbose('Reading Markers File ''%s''',10,markersfile)
    [markersfile_hasheader,~]=check_if_has_header(markersfile,3);

    % read markers data
    fid=fopen(markersfile);
    if fid < 0
        throw(MException('snp:make_D_from_seg:markersfileOpenError',...
                    'Error opening markers file ''%s''.',markersfile));
    end
    try
        markers = textscan(fid,'%q%q%q','ReturnOnError',0,'Headerlines',0+markersfile_hasheader);
    catch me
        fclose(fid);
        throw(MException('snp:make_D_from_seg:SegmentInputScanError',...
                        strcat('Error scanning markers file ''%s'':\n',me.message),...
                        markersfile));
    end
    fclose(fid);
    
    % translate marker chromosome column to numbers
    mrkr_chrn = double(chromosome2num(markers{2}));
    % if errors, generate message about the first error
    if any(isnan(mrkr_chrn))
        badline = find(isnan(mrkr_chrn),1,'first');
        throw(MException('snp:make_D_from_seg:BadMarkers',...
            'Invalid chromosome ''%s'' in markers file ''%s'',line %d',...
            markers{2}{badline},markersfile,badline+markersfile_hasheader));
    end
    % translate marker positions to numbers
    mrkr_pos = str2double(markers{3});
    % if errors, generate message about the first error
    if any(isnan(mrkr_pos))
        badline = find(isnan(mrkr_pos),1,'first');
        throw(MException('snp:make_D_from_seg:BadMarkers',...
            'Non-numeric position ''%s'' in markers file ''%s'', line %d',...
            markers{3}(badline),markersfile,badline+markersfile_hasheader));
    end
else
    % markers "file" is a cell array of columns
    verbose('Using supplied markers',20);
    if length(markersfile)==2 % no marker names
        tmp={'Empty'};
        markers{1}=tmp(ones(length(markersfile{1}),1));
        %  cellstr(num2str((1:length(markersfile{1}))'));
        markers(2:3)=markersfile(1:2);
    else
        markers=markersfile(1:3);
    end
    % get chromosome and base position
    if iscellstr(markers{2})
        mrkr_chrn = double(chromosome2num(markers{2}));
    else
        mrkr_chrn = markers{2};
    end
    if iscellstr(markers{3})
        mrkr_pos = str2double(markers{3});
    else
        mrkr_pos = markers{3};
    end
end

% sort markers by genomic position
markerpos=double(mrkr_chrn)*1e11+double(mrkr_pos);
if ~issorted(markerpos)
    verbose('Markers in markersfile require sorting!');
    [markerpos,order] = sort(markerpos);
    mrkr_chrn = mrkr_chrn(order);
    mrkr_pos = mrkr_pos(order);
end
%! TODO combine with next operation to make more efficient 

% remove duplicate markers
[upos,ui] = lunique(markerpos,'first');
if length(markerpos)~=length(upos)
    verbose('Non-unique positions ... using first marker from each position',10);
end

% choice of how to process markers...
if options.compress_markers
    % ui will compress markers; segmap will not remap segments
    segmap = 1:length(upos);
else
    % segmap will remap segments; ui will not compress markers
    segmap = ui;
    ui = 1:length(markerpos);
end

%% read segmented data file
[segs,hasheader] = read_segfile(segfile);

% if not specified, auto detect units (log ratio or copy number)
detect_logR = min(segs.copyN) < 0;
if ~isfield(options,'islog') || isempty(options.islog)
    options.islog = detect_logR;
elseif detect_logR == ~options.islog
    if detect_logr
        warning('data appears to be log ratio, option species copy units');
    else
        warning('data appears to be copy units, user option specifies log units');
    end
end

% create sample index, match all samples to unique samples
sample_names = unique_keepord(deblank(segs.sample),'stable');

% create sample index, match all samples to unique samples
[~,samplex] = match_string_sets_hash(sample_names,segs.sample);
samplex = samplex(:); % (for single sample case)


% find marker(s) for segment start positions
[matchbeg,st1] = match_num_sets(upos,segs.chrn*1e11+double(segs.start));
mismatch_beg = full(~any(matchbeg));
% find marker(s) for segment end positions
[matchend,en1] = match_num_sets(upos,segs.chrn*1e11+double(segs.end));
mismatch_end = full(~any(matchend));

if any(mismatch_beg) || any(mismatch_end)
    Nmismatch = sum(mismatch_beg) + sum(mismatch_end);
    badline = find(mismatch_beg|mismatch_end,1,'first');
    if mismatch_beg(badline)
        badpos = segs.start(badline);
    else
        badpos = segs.end(badline);
    end
    throw(MException('snp:make_D_from_seg:BadSegment',...
        ['%d segment start or end positions in ''%s'' do not match any markers in ''%s''.\n',...
        'First bad position is %s:%d at line %d.'],...
        Nmismatch,segfile,markersfile,segs.chr{badline},badpos,badline+hasheader));
end


%% check for overlapping segments

% get ordering vector for sort by (1) sample (2) chromosome (3) segment start
[~,order] = sort(st1);
[~,neord] = sort(segs.chrn(order));
order = order(neord);
[~,neord] = sort(samplex(order));
order = order(neord);

% order each segment's begin and end positions
begends = [st1(order) en1(order)];
[begends,colord] = sort(begends,2);
flipped = colord(:,1)==2;
if any(flipped)
    % warn about reversed segments
    firstflip = find(flipped,1,'first') + hasheader;
    warning( 'snp:make_D_from_seg:ReverseSegment',...
            ['%d segments are reversed in file ''%s''.\n' ...
             'The first reversed segment is on line %d.'], ...
            sum(flipped), segfile, firstflip);
end

% figure out which segment ends are interior to chromosome/samples 
inner_segends = (diff(samplex(order))==0)&(diff(segs.chrn(order))==0);
% get the overlap sizes within sample/chromosome (negative = gap)
osizes = inner_segends .* (begends(1:end-1,2) - begends(2:end,1) + 1);

% compare segment end positions with begin positions of next segment
if any(osizes > 0)
    if max(osizes) > 1
        % error if overlap sizes greater than one
        badseg = find(osizes > 1,1,'first');
        badlines = order([badseg,badseg+1]) + hasheader;
        throw(MException('snp:make_D_from_seg:segOverlap',...
                    ['%d segment overlaps detected in file ''%s''.\n',...
                    'First overlap detected between segments at lines %d and %d.'],...
                    sum(osizes>1),segfile,badlines(1),badlines(2)));
    else
        % tolerate single probe overlaps by shortening segment ends
        olap1s = [osizes == 1;false];
        warning('snp:make_D_from_seg:zeroLengthSegment',...
                'Shortened %d segments in ''%s'' that overlap by one marker.',...
                sum(olap1s),segfile);
        begends(olap1s,2) = begends(olap1s,2) - 1;
        % error if any segments are shortened to zero length
        if any(begends(:,2) - begends(:,1) < 0)
            % find the original file lines with zero-shortened segments 
            badlines = sort(order(begends(:,2)-begends(:,1) < 0));
            warning('snp:make_D_from_seg:zeroLengthSegment',...
                    ['Removing %d segments shortened to length zero in file ''%s''.\n' ...
                     'First segment shorted out is at line %d.'],...
                    length(badlines),segfile,badlines(1) + hasheader);
            % remove short-out segments
            st1(badlines) = [];
            en1(badlines) = [];
            samplex(badlines) = [];
            segs.copyN(badlines) = [];
        end
    end
end

%% validate copy number values

if all(isnan(segs.copyN))
    throw(MException('snp:make_D_from_seg:noCopyData',...
                     'No copy number data detected.'))
end

%% build the output data structure
D = struct;
D.sdesc = sample_names;
D.chrn = mrkr_chrn(ui);
D.pos = mrkr_pos(ui);
if options.marker_field
    D.marker = markers{1}(ui);
end
if options.chr_field
    D.chr = num2chromosome(D.chrn);
end
D.suppress_history = options.suppress_history;
D.islog = options.islog;
D.isMB = options.isMB;

% create data array
if options.use_segarray
    % create a SegArray for D.dat
    D.dat = SegArray.fromSegments(segmap(st1),segmap(en1),samplex,segs.copyN,...
                                  NaN,length(D.pos),length(D.sdesc));
else
    % create a full array for D.dat
    D.dat=nan(length(D.pos),length(D.sdesc));
    for i=1:length(D.sdesc)
        segidx=find(samplex==i); % all segments for sample i
        for j=1:length(segidx)
            s=segidx(j);
            D.dat(segmap(st1(s)):segmap(en1(s)),i)=segs.copyN(s);
        end
        if mod(i,50)==0
            disp(i);
        end
    end
end

%% subfunctions

