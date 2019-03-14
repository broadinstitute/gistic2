function D = make_D_from_segseq_data(segfile,options)
%MAKE_D_FROM_SEGSEQ_DATA make "D" data structure from markerless segmented CN data
%
%   D = make_D_from_segseq_data(SEGFILE,OPTIONS) 
%
% SEGFILE is a path to a file (string) containing segmented data
%    in whitespace delimited columns in this order: 
%        sample, chromosome, start base,end base, num markers, copy number
% OPTIONS is a structure containing optional values
%   OPTIONS.spacing is the maximum spacing between marker positions, 
%   default 10K base (must be between 1K base and 100K base)
%   OPTIONS.islog indicates that the data are log ratio (default
%   detect from data)
%

%% process input parameters

%{
if ~exist('spacing','var')
    spacing = 1e4;
end
if spacing < 1e3 || spacing > 1e5
    throw(MException('snp:make_D_from_segseq_data:bad_parameter',...
       'Maximum marker spacing parameter must be between 10kbase and 100kbase.'));
end
%}

% options
if ~exist('options','var')
    options = struct;
end

% provide defaults for undefined options
%! options = impose_default_value(options,'islog',true); %!!! allow detection
options = impose_default_value(options,'spacing',10000);

%% read segmented data file
[segs,hasheader] = read_segfile(segfile);

%% if not set as an option, detect if data is log scaled by presence of negative numbers
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

[mrkr_pos,mrkr_chrn] = space_markers(segs,options.spacing);
upos = double(mrkr_chrn)*1e11+double(mrkr_pos);

% create sample index, match all samples to unique samples
sample_names = unique_keepord(deblank(segs.sample));
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
D.chrn = mrkr_chrn;
D.pos = mrkr_pos;
D.suppress_history = true;
D.islog = options.islog;
D.isMB = false;

% create copy number data array
D.dat = SegArray.fromSegments(st1,en1,samplex,segs.copyN,NaN,length(D.pos),length(D.sdesc));

%% SUBFUNCTIONS

% interpolate markers to cover segments
function [mrkr_pos,mrkr_chrn] = space_markers(segs,spacing);
uchrn = unique(segs.chrn);
Nchrn = length(uchrn);

c_pos = cell(Nchrn,1);
c_chrn = cell(Nchrn,1);

for i = 1:Nchrn
    c = uchrn(i);
    onchr = segs.chrn==c;
    segbounds = unique([segs.start(onchr);segs.end(onchr)]);
    segaps = diff(segbounds);
    nchrmarks = sum(ceil(segaps/spacing))+1;
    chrmarks = nan(nchrmarks,1);

    m = 1;
    for j=1:length(segaps)
        gap = segaps(j);
        n = ceil(gap/spacing); % number of interpolated positions
        nspace = gap/n;

        chrmarks(m:m+n-1) = segbounds(j)+round((0:n-1)*nspace);
        m = m + n;
    end
    chrmarks(end) = segbounds(end);
    c_pos{c} = chrmarks;
    c_chrn{c} = repmat(c,nchrmarks,1);
end

mrkr_pos = cat(1,c_pos{:});
mrkr_chrn = cat(1,c_chrn{:});
