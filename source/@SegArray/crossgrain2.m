function [Q R] = crossgrain2(S,ccfh)
%CROSSGRAIN compression op with 2 return values over dim 2 of a SegArray 
%   [Q R] = crossgrain2(S,CCFH)
%
%   S is the input SegArray matrix
%   CCFH is a handle to the compression function, e.g. @(x) = sum(x,2);
%   Q R are the resulting vectors (or scalar if S is a rwo vector)
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

[M N] = size(S.bpts);
[K,J] = find(S.bpts);
if M == 1
    % return a scalar for row vectors
    [Q R] = ccfh(S.vals');
    R = J(R);
else
    % resulting breakpoints will be union of column breakpoints
    brkrows = unique(K);
    nrows = length(brkrows);
    % create output SegArrays
    % (know breakpoints, but values still need to be calculated)
    Q = SegArray();
    Q.bpts = sparse(brkrows,ones(nrows,1),1:nrows,M,1);
    Q.vals = nan(nrows,1);
    R = Q;
    % scale chunks according to SegArray memory setting
    chunksize = ceil(SegArray.MAX_INDICES / N);
    % initialize chunk per-row pre-indices
    prerowx = S.bpts(1,:) - 1;
    % loop over chunks
    for chunkstart = 1:chunksize:nrows
        chunkend = min(chunkstart+chunksize-1,nrows);
        chunk = double(logical(S.bpts(brkrows(chunkstart:chunkend),:)));
        chunk(1,:) = chunk(1,:) + prerowx;
        % convert breakpoints to value indices 
        chunk = cumsum(full(chunk),1);
        % update chunk per-row pre-indices
        prerowx = chunk(end,:);
        % perform cross compression operation on a chunk
        if size(chunk,1) == 1
            [Q.vals(chunkstart:chunkend) R.vals(chunkstart:chunkend)] = ccfh(S.vals(chunk)');
        else
            [Q.vals(chunkstart:chunkend) R.vals(chunkstart:chunkend)] = ccfh(S.vals(chunk));
        end
    end
end
