function verbosedisp(var,thresh)
%VERBOSEDISP  Display "whos" of input variables if VERBOSE_LEVEL is set to
%thresh.
%
% verbosedisp(VAR)
%
%   See also:  SET_VERBOSE_LEVEL, VERBOSE

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.



global VERBOSE_LEVEL

if VERBOSE_LEVEL >= thresh

    disp([inputname(1) ' = '])
    disp(var)
    
end
