function S = cat(DIM,varargin)
% CAT concatenate SegArrays with full arrays or SegArrays
% The result is always a segarray

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

switch DIM
  case 1
    S = vertcat(varargin{:});
  case 2
    S = horzcat(varargin{:});
  otherwise
    throwAsCaller(MException('MATLAB:SEGARRAY:dimTooHighMan', ...
               'SegArray does not catenate along dimensions higher than 2.'));
end

% TODO: (1) test
% TODO: (2) try...catch...throw as caller after horzcat and vertcat are tested enough
