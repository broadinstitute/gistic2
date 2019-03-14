function S = fromSegments(B,E,C,V,F,M,N,TOL)
% SEGARRAY.FROMSEGMENTS create a SegArray from a list of segments.
%
%   S = SegArray.FROMSEGMENTS(B,E,C,V,F,M,N) creates an M-by-N SegArray
%   using vectors of begin and end rows (B and E) for each segment
%   along with the corresponding columns in C and values in V. Any
%   gaps between segments are filled with the scalar value F. 
%  
%   Segments that extend beyond the boundaries of an M-by-N matrix
%   will generate a MATLAB:SEGARRAY:segOvershoot exception. The
%   order of B and E is unimportant.
%   
%   The B and E row subscript arguments are mandatory and must be vectors
%   of equal length; the remaining argumentare optional with default values:
%       
%   C can be a vector or scalar value to repeat with a default of
%   one (creating a column segvector).
%
%   V can be a vector of values or repeated scalar, default TRUE.
%
%   F must be a scalar filler value. If V is logical, the filler default 
%   is FALSE, otherwise the default is NaN.
%
%   If M and/or N are omitted, then a minimum size SegArray is created
%   that fits the specified segments: M defaults to MAX([B;E]) and N
%   defaults to max(C).
%
%   If the test overlap argument TOL is present and nonzero, then
%   the overlap of any input segments will cause a
%   MATLAB:SEGARRAY:segOverlap exception to be thrown. The default
%   behavior is for segments defined later in the lists to
%   overwrite earlier ones. 

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% validate mandatory segment row begin/end arguments
if ~exist('B','var') || isempty(B) || ~exist('E','var') || isempty(E)
    throwAsCaller(MException('MATLAB:SEGARRAY:missingArg',...
                             'SEGARRAY.fromSegment requires begin and end row arguments.'));
end
B = B(:);
E = E(:);
nsegs = length(B);
if length(E) ~= nsegs
    throwAsCaller(MException('MATLAB:SEGARRAY:argSizeMismatch',...
                             'SEGARRAY must have same number of begin and end rows.'));                    
end
% validate row index B
try
    validateattributes(B,{'numeric'},{'positive','integer'});
catch
    throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                             'SEGARRAY.fromSegment begin rows must be real positive integers.'));                    
end

%validate row index E
try
    validateattributes(E,{'numeric'},{'positive','integer'});
catch
    throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                             'SEGARRAY.fromSegment end rows must be real positive integers.'));                    
end

% validate columns in C
if ~exist('C','var') || isempty(C)
    C = 1;
    N = 1;
else
    C = C(:);
    % validate column index C
    try
        validateattributes(C,{'numeric'},{'positive','integer'});
    catch
        throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                                 'SEGARRAY.fromSegment columns must be real positive integers.'));                    
    end
    if ~isscalar(C)
        if length(C) ~= nsegs
            throwAsCaller(MException('MATLAB:SEGARRAY:argSizeMismatch',...
                                     'SEGARRAY column vector size mismatch.'));                    
        end
    end
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
% validate filler value F
if ~exist('F','var') || isempty(F)
    if islogical(V)
        F = false;
    elseif ischar(V)
        F = ' ';
    else
        F = NaN;
    end
else
    if ~isscalar(F)
        throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                                 'SEGARRAY fill value must be a scalar.'));
    end
end
% validate array dimensions 
if ~exist('M','var') || isempty(M)
    M = max([B;E]);
else
    if ~isscalar(M)
        throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                                 'SEGARRAY dimension must be a scalar.'));
    end
    if any([B;E] > M)
        throwAsCaller(MException('MATLAB:SEGARRAY:segOvershoot',...
                                 'SEGARRAY row outside of array.'));
    end
end
if ~exist('N','var') || isempty(N)
    N = max(C);
else
    if ~isscalar(M)
        throwAsCaller(MException('MATLAB:SEGARRAY:invalidArg',...
                                 'SEGARRAY dimension must be a scalar.'));
    end
    if any(C > N)
        throwAsCaller(MException('MATLAB:SEGARRAY:segOvershoot',...
                                 'SEGARRAY column outside of array.'));
    end
end
% repeat anything that needs it
if isscalar(V)
    V = repmat(V,nsegs,1);
end
if isscalar(C)
    C = repmat(C,nsegs,1);
end
% default test overlap argument
if ~exist('TOL','var') || isempty(TOL)
    TOL = false;
end

% sort segments and polarize them (so B <= E)
[sx order] = sort( sub2ind([M N],B,C) );
V = V(order);
ex = sub2ind([M N],E,C);
segs = sort([sx ex(order)],2);
% overlap testing
if TOL
    if (nsegs > 1)
        if ~all(segs(2:end,1) > segs(1:end-1,2))
            throwAsCaller(MException('MATLAB:SEGARRAY:segOverlap',...
                                     'SEGARRAY segments overlap.'));
        end
    end
else
    %Clean up possibly overlapping data:
    % (1) trim ends of segments to make them strictly non-overlapped
    olaps = segs(1:end-1,2) >= segs(2:end,1);
    segs(olaps,2) = segs(1+find(olaps),1) - 1;
    % (2) eliminate nil segments and strictly overlapped segments
    badsegs = (segs(:,2) < segs(:,1)) | [false;segs(1:end-1,2) > segs(2:end,2)];
    segs(badsegs,:) = [];
    V(badsegs) = [];
    % (no need to compress B,E,C which are no longer used)
end

% co-sort starts, ends with column breaks
% priority (1) segment start values, (2) segment end NaNs, (3) column break NaNs
values = [V; repmat(F,nsegs,1); repmat(F,N,1)];
idx = [segs(:,1); ...
       segs(:,2) + 1; ...
       sub2ind([M N],ones(N,1),(1:N)')];
[idx order] = sort(idx);
% line up values and mark the first ones
values = values(order);
undup = ne(idx, circshift(idx,1)) & idx <= M*N;
undup(1) = true;
% put breakpoints and their values in a SegArray
S = SegArray;
S.vals = values(undup);
[i j] = ind2sub([M N], idx(undup));
S.bpts = sparse(i,j,(1:size(S.vals))',M,N); % LOGBPTS
%fin
