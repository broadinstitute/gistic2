%% GET_NUMVALS returns the number of segment breakpoints stored in a SEGARRAY, and will be renamed soon

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function num_vals = get_numvals(S)
    num_vals = nnz(S.bpts);
end
