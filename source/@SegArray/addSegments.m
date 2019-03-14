function S = addSegments(S,B,E,C,V)
% SEGARRAY.ADDSEGMENTS adds a list of segments to a numeric SegArray.
%
%   S = SegArray.ADDSEGMENTS(S,B,E,C,V) adds segments to a SegArray S
%   using vectors of begin and end rows (B and E) for each segment
%   along with the corresponding columns in C and values in V. 
%  
%   Segments that extend beyond the boundaries of S will generate a
%   MATLAB:SEGARRAY:segOvershoot exception. The order of B and E is
%   unimportant. 
%   
%   The B and E row subscript arguments are mandatory and must be vectors
%   of equal length; the remaining argumentare optional with default values:
%       
%   C can be a vector or scalar value to repeat with a default of
%   one (creating a column segvector).
%
%   V can be a vector of values or repeated scalar, default TRUE.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% validate mandatory segment row begin/end arguments
if ~exist('B','var') || ~exist('E','var')
    throwAsCaller(MException('MATLAB:SEGARRAY:missingArg',...
                             'SegArray.addSegments requires begin and end row arguments.'));
end
B = B(:);
E = E(:);
nsegs = length(B);

% get SegArray size
[M N] = size(S.bpts);

if length(E) ~= nsegs
    throwAsCaller(MException('MATLAB:SEGARRAY:argSizeMismatch',...
                             'SEGARRAY must have same number of begin and end rows.'));
end

% allow empty segments (null operation) 
if nsegs == 0
    return;
end

% validate row index B
try
    validateattributes(B,{'numeric'},{'positive','integer'});
catch me
    throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                             'SegArray.addSegments begin rows must be real positive integers.'));
end

%validate row index E
try
    validateattributes(E,{'numeric'},{'positive','integer'});
catch me
    throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                             'SegArray.addSegments end rows must be real positive integers.'));
end

% validate columns in C
if ~exist('C','var') || isempty(C)
    C = 1;
else
    C = C(:);
    % validate column index C
    try
        validateattributes(C,{'numeric'},{'positive','integer'});
    catch me
        throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                                 'SegArray.addSegments columns must be real positive integers.'));
    end
    if ~isscalar(C)
        if length(C) ~= nsegs
            throwAsCaller(MException('MATLAB:SEGARRAY:argSizeMismatch',...
                                     'SEGARRAY column vector size mismatch.'));
        end
    end
end

% validate segment ranges vs SegArray size
if any([B;E] > M)
    throwAsCaller(MException('MATLAB:SEGARRAY:segOvershoot',...
                             'SEGARRAY row outside of array.'));
end

if any(C > N)
    throwAsCaller(MException('MATLAB:SEGARRAY:segOvershoot',...
                             'SEGARRAY column outside of array.'));
end

% validate segment values in V
if ~exist('V','var') || isempty(V)
    V = true;
else
    if ~isscalar(V)
        V = V(:);
        if length(V) ~= nsegs
            throwAsCaller(MException('MATLAB:SEGARRAY:argSizeMismatch',...
                                     'SEGARRAY value vector size mismatch.'));                    
        end
    end
end

% repeat anything that needs it
if isscalar(V)
    V = repmat(V,nsegs,1);
end
if isscalar(C)
    C = repmat(C,nsegs,1);
end

% polarize BE segments
sx = sub2ind([M N],B,C);
ex = sub2ind([M N],E,C);
segs = sort([sx ex],2);

% combine new segment values with existing segment values
combivals = [V; S.vals];
% Collect possibly overlapping new breakpoints from combined segments
% Exactly twice the size of combined values: all segments starts, followed 
% by all segment ends. The last breakpoint is a pseudo break just over the
% boundary of the SegArray that is used to extend the last existing
% segment.
Sbrks = find(S.bpts);
newbpts = [segs(:,1); Sbrks; segs(:,2)+1; Sbrks(2:end); M*N+1];

% sort new breakpoints, overlappers will be grouped
[newbpts order] = sort(newbpts);
% mark beginning of groups of identical indices
starts = [true;0~=diff(newbpts)];
newbpts = newbpts(starts);
% create N:1 mapping from combined value position to new value position
bpidx = cumsum(starts);
backmap = zeros(size(order));
backmap(order) = bpidx;
% create output object
S = SegArray;

% accumulate segment values in S.vals
nsgs = length(combivals);
%! use Mex routine to sum ranges fast
S.vals = SegArray.sumranges(backmap(1:nsgs),backmap(nsgs+1:end)-1,combivals,bpidx(end));
%! the code above is equivalent to the commented out code below
%{
S.vals = zeros(bpidx(end),1);
for k=1:length(combivals)
    vidx = backmap(k):backmap(k+length(combivals))-1;
    S.vals(vidx) = S.vals(vidx) + combivals(k); %!!! > 99% of time spent here
end
%}

% remove final pseudo-breakpoint and its value 
newbpts(end)=[];
S.vals(end) = [];
% make 2D sparse array for breakpoints from linear indices
[i j] = ind2sub([M N], newbpts);
S.bpts = sparse(i,j,(1:size(S.vals))',M,N); % LOGBPTS
%fin
