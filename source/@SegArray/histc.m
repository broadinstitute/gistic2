function [N,BIN] = histc(S,EDGES,DIM)
% SegArray implementation of HISTC 
%   Returned N is a full vector/array
%   Returned BIN is a SegArray

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


    if ~exist('DIM', 'var')
        DIM = 1;
    end
    if DIM == 1 && exist('EDGES','var') && isvector(EDGES) && ...
            issorted(EDGES) && ~isempty(EDGES) && isa(S,'SegArray')
        % sort the segment values
        [sorted,order] = sort(S.vals);
        % calculate segment sizes
        idx = find(S.bpts);
        sizes = [idx(2:end);numel(S.bpts)+1]-idx;
        % line them up
        sizes = sizes(order);
        [~,J] = find(S.bpts);
        J = J(order);
        % nested binning loop
        I = zeros(size(sorted));
        bin = 1;
        k = find(sorted >= EDGES(1),1);
        if ~isempty(k)
            while k <= length(I) && sorted(k) <= EDGES(end)
                while  bin <= numel(EDGES)-1 && sorted(k) >=  EDGES(bin+1)
                    bin = bin + 1;
                end
                I(k) = bin;
                k = k + 1;
            end
        end
        mask = logical(I);
        N = accumarray([I(mask) J(mask)],sizes(mask),[length(EDGES) size(S.bpts,2)]);
        % 2nd output is SegArray with values replaced by bin indices
        if nargout > 1
            BIN = S;
            BIN.vals(order) = I;
        end
    else
        % convert to full for other dimensions
        try
            if nargout > 1
                [N,BIN] = histc(full(S),full(EDGES),DIM);
            else
                N = histc(full(S),full(EDGES),DIM);
            end
        catch me
            throwAsCaller(me);
        end
    end
