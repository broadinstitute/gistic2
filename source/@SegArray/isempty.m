% SegArray implementation of ISEMPTY.
%    TF = ISEMPTY(A)
% Emulates normal behavior of ISEMPTY.
%    (see HELP ISEMPTY for details)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function TF = isempty(A)
    try
        TF = isempty(A.bpts); % queries shape
        if TF
            segwarn(A,'Empty SegArray allowed to exist!');
        end 
    catch me
        throwAsCaller(me);
    end
end
