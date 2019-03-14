function R = crossgrain(S,ccfh)
%CROSSGRAIN execute a compression operation over dimension 2 of a SegArray
%
%   R = crossgrain(S,CCFH)
%
%   S is the input SegArray matrix
%   CCFH is a handle to the compression function, e.g. @(x) = sum(x,2);
%   R is the resulting vector (or scalar if S is a rwo vector)
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

[M N] = size(S.bpts);

if M == 1
    % return a scalar for row vectors
    R = ccfh(S.vals');
else
    % resulting breakpoints will be union of column breakpoints
    [i,~] = find(S.bpts);
    brkrows = unique(i);
    nrows = length(brkrows);
    % create output SegArray
    % (know breakpoints, but values still need to be calculated)
    R = SegArray();
    R.bpts = sparse(brkrows,ones(nrows,1),1:nrows,M,1);

    % set the type of the R.vals according to function results 
    if islogical(ccfh(full(S.vals(1,:))))
        R.vals = false(nrows,1);
    else
        R.vals = nan(nrows,1);
    end
    
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
            R.vals(chunkstart:chunkend) = ccfh(S.vals(chunk)');
        else
            R.vals(chunkstart:chunkend) = ccfh(S.vals(chunk));
        end
    end
end
