function B = double(A)
%SegArray implementation of DOUBLE.
%    B = DOUBLE(A)
% If A is a SegArray of double values, assume a caller is trying
% to cast SegArray into something it understands (e.g. subsasgn).
% Otherwise, just cast type of SegArray to double. (Sorry this is 
% so weird.)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if strcmp(class(A.vals),'double')
    B = full(A);
else
    B = A;
    try
        B.vals = double(A.vals);
    catch me
        throwAsCaller(me);
    end
end
