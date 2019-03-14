%SUM - overloaded sum for SegArray
%
%   R = SUM(S)
%   R = SUM(S,DIM)
%
% Return value is a full row vector if DIM = 1, a SegArray column
% vector if DIM = 2

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% Note: currently maintained in parallel with nansum.m
function R = sum(S,dim)
    if ~exist('dim','var') || isempty('dim')
        dim = 1 + (1 == size(S.bpts,1));
    else
        try
            % SegArray does not support dimensions higher than two
            validateattributes(dim,{'numeric'},{'positive','integer','scalar'});
        catch
            % generate error for invalid dim values
            try
                sum(S.bpts,dim);
                errNoError(S);
            catch me
                throwAsCaller(me);
            end
        end
    end

    switch dim
      case 1 %% sum with the grain of SegArray
        if isempty(S)
            R = zeros(1,size(S,2));
        else
            % The idea for summing "with the grain" is to do the calculation
            % with values and counts of those values.
            [M,N] = size(S.bpts);
            if M == 1
                R = full(S);
            else
                [I,J] = find(S.bpts);
                L = sub2ind([M N],I,J);
                counts = [diff(L);numel(S.bpts)-L(end)+1];
                wgtvals = S.vals .* counts;
                R = full(sum(sparse(I,J,wgtvals,M,N)));
            end
        end
      case 2 %% sum across the grain of SegArray
        if isempty(S)
            R = zeros(size(S,1),1);
        else
            if size(S.bpts,1) == 1
                R = sum(S.vals);
            else
                R = crossgrain(S,@(x) sum(x,2));
            end
        end
      otherwise % dim neither 1 nor 2
        throwAsCaller(MException('MATLAB:SEGARRAY:dimTooHighMan', ...
                  'Indexing higher than two dimensions unsupported by SegArray objects'));
    end
end
