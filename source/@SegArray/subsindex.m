% Overloaded subsindex function for SegArray objects

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function I = subsindex(A)
    if islogical(A.vals)
        I = find(full(A))-1;
    else
        I = full(A) - 1; 
    end
    % (can't do much else without access to the indexed object...)
end
