function R = any(S,dim)
%ANY - overloaded ANY (OR summary) for SegArray
%
%   R = ANY(S)
%   R = ANY(S,DIM)
%
% Return value is a full row vector if DIM = 1, 
% a SegArray column vector if DIM = 2 and S is 
% not a row vector.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


    [M N] = size(S.bpts);
    
    if ~exist('dim','var') || isempty(dim)
        if M || N
            dim = find(size(S)~=1,1);
            % if scalar operand, return non-SegArray result
            if isempty(dim)
                R = logical(full(S));
                return
            end
        else
            % any([]) is false
            R = false;
            return
        end
    end

    if dim == 1 %% compress across columns ("with the grain")
        if isempty(S.bpts)
            R = zeros(1,N);
        else
            % optimize row and column vectors
            if M == 1
                R = logical(S.vals)';
            elseif N == 1
                R = any(S.vals);
            else
                % Testing "with the grain" for any nonzeros
                [~,J] = find(S.bpts);
                R = accumarray(J(:),logical(S.vals),[N 1],@any,false,false)';
            end
        end
    
    elseif dim == 2 %% compress across rows ("against the grain")
        if isempty(S.bpts)
            R = zeros(M,1);
        else
            R = crossgrain(S,@(x) any(x,2));
        end
    else
        R = logical(S);
    end

