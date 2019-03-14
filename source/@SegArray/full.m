% FULL Convert SEGARRAY compressed segmented array to full array

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function M = full(S)
    [m n] = size(S.bpts);
    nzs = reshape(logical(S.bpts),numel(S.bpts),1); % LOGBPTS
    ii = reshape(cumsum(full(nzs)),m,n);
    M = S.vals(ii);
    % adapt to Matlab's special behavior for row vector indexing
    if m == 1
        M = M';
    end
