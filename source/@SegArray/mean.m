function y = mean(x,dim)
% OVERLOADED mean function for SegArray object.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if nargin == 1
    dim = min(find(size(x) ~=1));
end

if isempty(dim)
    dim = 1;
end

y = sum(x,dim) ./ size(x,dim);
