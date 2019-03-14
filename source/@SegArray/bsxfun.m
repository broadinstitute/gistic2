function C = bsxfun(fun,A,B)
% SegArray bsxfun helper

% GISTIC software version 2.0
% Copyright (c) 2011-2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steven Schumacher, Nicolas Stransky, 
% Mike Lawrence, Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

erred = false; % innocent until proven guilty...

if ~exist('fun','var') || isempty(fun) || ~isa(fun,'function_handle')
    error('first arg must be a function handle');
end

if isscalar(A) || isscalar(B) || all(size(A) == size(B))
    C = binary(A,B,fun);
elseif size(A,1) == size(B,1)
    if size(A,2) == 1
        A = repmat(A,1,size(B,2));
    elseif size(B,2) == 1;
        B = repmat(B,1,size(A,2));
    else
        % signal error
        erred = true;
    end
    C = binary(A,B,fun);
elseif size(A,2) == size(B,2)
    if size(A,1) == 1
        A = repmat(SegArray(A),size(B,1),1);
    elseif size(B,1) == 1;
        B = repmat(SegArray(B),size(A,1),1);
    else
        % signal error
        erred = true;
    end
    C = binary(A,B,fun);
else
    % signal error
    erred = true;
end

if erred
    % throw artificially generated non-singleton dimansion mismatch error
    try
        bsxfun(@plus,magic(3),1:4)
    catch me
        throwAsCaller(me);
    end
end


