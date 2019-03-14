function res = isunique(x,do_rows)
%ISUNIQUE Test for uniqueness of elements (or rows) in a vector (or matrix)
%
% RES = isunique(X,DO_ROWS)
%
% returns TRUE if every element of X is unique. If DO_ROWS is provided and
% true, then rows of matrix X are tested for uniqueness.

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if size(x,1)==1
  x = x';
end

if exist('do_rows','var') & do_rows
  u = unique(x,'rows');
else
  u = unique(x);
end

res = (size(u,1)==size(x,1));
