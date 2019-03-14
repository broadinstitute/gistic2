% internal SEGARRAY warning message function
% issues provided warning message if global SEGWARNING is defined
% S is a dummy argument to make it a class method

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function segwarn(S,msg)
    global SEGWARNING
    if ~isempty(whos('global','SEGWARNING')) & SEGWARNING
        warning(['(SEGARRAY) ' msg]);
    end
end
