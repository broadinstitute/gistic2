function [I,J,V] = find(X,K,str)
% Overloaded FIND function for SegArray objects.
% Always returns a non-Segarray

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  
  if ~exist('K','var') || isempty(K)
    K = numel(X);
  elseif ~isnumeric(K)
    error('K must be a numeric value');
  end
  
  if ~exist('str','var') || isempty(str)
      str = 'first';
  elseif ~strcmp(str,'first') && ~strcmp(str,'last')
      error('Must input either "first" or "last" for str');
  end
  
  vals = X.vals;
  bpts = X.bpts;

  idx = find(bpts);

  % allocate array to store indices in
  nonzeros = sum(diff([idx(:); numel(X.bpts)+1]) .* (X.vals ~=0));
  I = zeros(nonzeros,1);
  % loop over segments
  fidx = 1;
  for j=1:length(idx)-1
    % store index range for each nonzero segment
    if vals(bpts(idx(j)))
      segsize = idx(j+1) - idx(j) ;
      I(fidx:fidx+segsize-1) = idx(j):idx(j+1)-1;
      fidx = fidx + segsize;
    end
  end
  if vals(bpts(idx(end)))
    I(fidx:nonzeros) = idx(end):numel(bpts);
  end
  %  display(size(I));%!!! test
  if K < length(I)
    switch str
     case 'first'
      I = I(1:K);
     case 'last'
      I = I(end-K+1:end);
    end
  end
  
  if nargout > 3
      error('Only 3 output variables supported...');
  end
  
  if size(X.bpts,1) == 1
      I = I';
  end
  
  if nargout == 3
      V = subsref(X,substruct('()',{I}));
  end

  if nargout >= 2
      II = I;
      I = mod(II,size(X,1));
      I(I == 0) = size(X,1);
      J = ceil(II/size(X,1));
  end
  
  
