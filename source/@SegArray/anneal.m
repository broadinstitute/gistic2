function S = anneal(S,supersegs)
%ANNEAL - Joins any equal-valued adjacent segments in a SegArray.
%   S = ANNEAL(S,SUPERSEGS)
% Compacts the storage of the SegArray object S by combining adjacent 
% segments that have the same value. Column header breakpoints and,
% optionally, row breakpoints specified in the vector SUPERSEGS, are preserved.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% keep non-identical values
keep = [true;diff(S.vals)~=0]; % was: (S.vals ~= circshift(S.vals,1));
nans = isnan(S.vals);
% only keep first NaN in a run
keep(nans & circshift(nans,1)) = false;
keep(1) = true;
% keep column breakpoints - COLBRK
[I,J] = find(S.bpts);
I = I(:); %! if S.bpts is a row vector, so are I,J !!!
J = J(:);
keep = keep | (I == 1);
% also keep any breakpoints along the supersegments
% (note that this preserves but does not create supersegments)
if exist('supersegs','var') && ~isempty(supersegs)
    supidx = repmat(supersegs(:),1,size(S.bpts,2)) + ...
             repmat(0:size(S.bpts,1):prod(size(S.bpts))-1, length(supersegs), 1);
    keep = keep | ismember(I,supidx);
end
% pitch the rest
S.vals = S.vals(keep);
I = I(keep);
J = J(keep);
S.bpts = sparse(I,J,1:size(S.vals),size(S.bpts,1),size(S.bpts,2));
