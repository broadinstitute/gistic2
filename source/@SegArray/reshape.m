% SEGARRAY implementation of RESHAPE

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function R = reshape(S,varargin)
    switch length(varargin)
      case 2
        [m n] = varargin{:};
        % deal with empty arguments
        if isempty(m)
            if isempty(n)
                throwAsCaller(MException('MATLAB:getReshapeDims:unknownDim', ...
                      'Size can only have one unknown dimension.'));
            end
            if 0 == mod(numel(S.bpts),n)
                m = numel(S.bpts) / n;
            else
                throwAsCaller(MException('MATLAB:getReshapeDims:notDivisible', ...
                      'Product of known dimensions, %d, not divisible by total number of elements, %d.', ...
                      n, numel(S.bpts) ));
            end
        elseif isempty(n)
            if 0 == mod(numel(S.bpts),m)
                n = numel(S.bpts) / m;
            else
                throwAsCaller(MException('MATLAB:getReshapeDims:notDivisible', ...
                      'Product of known dimensions, %d, not divisible by total number of elements, %d.', ...
                      m, numel(S.bpts) ));
            end
        elseif m * n ~= numel(S.bpts)
            throwAsCaller(MException('MATLAB:getReshapeDims:notSameNumel', ...
                      'To RESHAPE the number of elements must not change.'));
        end
        % call 2-D reshape
        R = reshape2D(S,m,n);
      case 0
        throwAsCaller(MException('MATLAB:minrhs', ...
                      'Not enough input arguments.'));
      case 1
        v = num2cell(varargin{1});
        if length(v) == 2
            [m n] = v{:};
            if m * n == numel(S.bpts)
                R = reshape2D(S,m,n);
            else
                throwAsCaller(MException('MATLAB:getReshapeDims:notSameNumel', ...
                          'To RESHAPE the number of elements must not change.'));
            end
        else
            throwAsCaller(MException('MATLAB:getReshapeDims:sizeVector', ...
                      'Size vector must have at least two elements.'));
        end
      otherwise
        throwAsCaller(MException('MATLAB:SegArray:exactlyTwoDims', ...
                      'SEGARRAY does not support more than two dimensions'));
    end
end

%% inner function - two dimensional reshape
function R = reshape2D(S,m,n)
    R = SegArray();
    [i j] = find(S.bpts);
    idx = sub2ind(size(S.bpts),i(:),j(:));
    % merge adjacent segments with identical values
    values = S.vals;
    mask = (values ~= circshift(values,1));
    mask(1) = true;
    idx = idx(mask);
    values = values(mask);
    % establish column header breakpoints - COLBRK
    ndx = sub2ind([m n],ones(n,1),(1:n)');     % interleave column headers with data indices
    idxndx = [idx;ndx];
    [sorted,order] = sort(idxndx);
    dmask = (sorted ~= circshift(sorted,1));   % unique indices mask
    dmask(1) = true;
    propidx = cumsum( order <= length(idx) );  % propagate data indices
    R.vals = values( propidx(dmask) );         % duplicate values
    [i j] = ind2sub([m n],(sorted(dmask)));    % reshape breakpoints
    R.bpts = sparse(i,j,1:size(R.vals),m,n);
end
