function set_verbose_level(level)
% SET_VERBOSE_LEVEL(LEVEL) sets the global VERBOSE_LEVEL  to LEVEL.
%
%    Standard levels are:
%            0 no messages
%           10 general messeges 
%           20 more details
%           30 all details
%
%    See also VERBOSE

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.



global VERBOSE_LEVEL

VERBOSE_LEVEL=level;
