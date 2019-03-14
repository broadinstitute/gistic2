function S = zeros(M,N);
% Create an all-zeros SegArray
%   Equivalent to ZEROS, except limited to two-dimensional arrays

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    if ~exist('M','var')
        M = 1;
    end
    if isempty(M)
        M = 0;
    end
    if ~exist('N','var') || isempty(N)
        if isscalar(M)
            N = M;
        else
            if size(M,1) == 1 && size(M,2) >= 2
                N = M(2);
                M = M(1);
            else
                % let matlab generate error
                try
                    zeros(M,N);
                catch me
                    throwAsCaller(me);
                end
                errNoError(SegArray);
            end
        end
    end
    if isempty(N)
        N = 0;
    end
    % make all-zeros SegArray
    S = SegArray;
    if M == 0 || N == 0
        % empty SegArray
        S.vals = [];
        S.bpts = sparse(M,N);
    else
        % non-empty SegArray
        S.vals = zeros(N,1);
        S.bpts = sparse(ones(1,N),1:N,1:N,M,N);
    end
