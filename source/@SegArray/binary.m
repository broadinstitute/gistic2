% internal helper function for SEGARRAY

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function C = binary(A,B,op)

if ~exist('op','var') || isempty(op)
    error('Must supply operator for binary function!');
end

if ~isscalar(A) & ~isscalar(B) & any(size(A) ~= size(B))
    % TODO throwAsCaller
    error('Matrix dimensions must agree.');
end

if ~isscalar(A) && ~isa(A,'SegArray')
    A = SegArray(A);
end

if ~isscalar(B) && ~isa(B,'SegArray')
    B = SegArray(B);
end

C = SegArray;
%C = set_super_seg_starts(C,union(A.super_seg_starts,B.super_seg_starts));
if isscalar(A) && isscalar(B)
    % should not happen for proper call from SEGARRAY method
    segwarn(SegArray,'scalar op scalar'); 
    C.vals = op(full(A),full(B));
    C.bpts = sparse(1);
elseif isscalar(A) && ~isscalar(B)
    C.bpts = B.bpts;
    C.vals = op(full(A),B.vals);
elseif isscalar(B)
    C.bpts = A.bpts;
    C.vals = op(A.vals,full(B));
else
    if isempty(A.bpts)
        % if operands are empty, result is just operand shape/type
        C.bpts = A.bpts;
        C.vals = A.vals + B.vals;
    else
        % non-empty: co-sort indices
        jA = find(A.bpts);
        jB = find(B.bpts);
        [sorted order] = sort([jA(:);jB(:)]);
        maskA = ( order <= length(jA) );
        jjA = cumsum(maskA);
        jjB = cumsum(~maskA);
        % remove others that are duplicates from each set
        uniq = sorted ~= circshift(sorted,1);
        uniq(1) = true;
        jjA = jjA(uniq);
        jjB = jjB(circshift(uniq,-1));
        % save resulting merged breakpoints
        [m n] = size(A.bpts);
        [i j] = ind2sub([m n], sorted(uniq));
        C.bpts = sparse(i,j,1:length(jjA),m,n);
        % perform operation on vectors of values
        C.vals = op( A.vals(jjA), B.vals(jjB) );
        C = anneal(C);
    end
end
