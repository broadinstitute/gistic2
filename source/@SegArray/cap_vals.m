function S = cap_vals(S,caps)
%CAP_VALS range-limit the values in SegArray S
%   S = CAP_VALS(S,CAPS) efficiently performs the equivalent of
%   S(S > caps(1)) = caps(1) and S(S < caps(2)) = caps(2)
%   If only maximum caps(1) is given, caps(2) is set to -caps(1)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('caps','var') || isempty(caps)
    throwAsCaller(MException('MATLAB:SEGARRAY:missingArg',...
        'Must supply cap values when using cap_vals function.'));
end
if length(caps) == 1
    caps = [caps -1*caps];
end
S.vals(S.vals > caps(1)) = caps(1);
S.vals(S.vals < caps(2)) = caps(2);
S = anneal(S);
