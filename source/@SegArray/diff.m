function Y = diff(X,n,dim)
% Overload DIFF function for SegArray objects.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~exist('n','var') || isempty(n)
    n = 1;
end

if ~exist('dim','var') || isempty(dim)
    dim = 1 + (size(X,1) == 1);
end

if n > 1
    Y = diff(diff(X,n-1,dim),1,dim);
else
    switch dim
      case 1
        % diff along the grain
        switch size(X,1)
          case 0
            Y = X;
          case 1
            Y = reshape([],0,size(X,2));
          otherwise
            new_vals = [1; diff(X.vals)];
            [i j] = find(X.bpts);
            yy = find(i-1);
            idx = (j(yy)-1)*(size(X,1)-1)+(i(yy)-1);
            Y = SegArray.zeros(size(X,1)-1,size(X,2));
            Y = subsasgn(Y,substruct('()',{idx}),new_vals(yy));
        end
      case 2
        % diff across the grain
        switch size(X,2)
          case 0
            Y = X;
          case 1
            Y = reshape([],size(X,1),0);
          otherwise
            Y = subsref(X,substruct('()', {':',2:size(X,2)} )) ...
                - subsref(X,substruct('()', {':',1:size(X,2)-1} ));
            % i.e. Y =  X(:,2:end) - X(:,1:end-1);
        end
      otherwise
        throwAsCaller(MException('MATLAB:SEGARRAY:dimTooHighMan', ...
                  'Indexing higher than two dimensions unsupported by SegArray objects'));        
    end
end
        
    
