function meds = nanmedian(S,dim)
%NANMEDIAN SegArray implementation of nanmedian function
%
%   The returned median vector is a row vector.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


sz = size(S);

% default dimension
if ~exist('dim','var') || isempty(dim)
    if sz(1) > 1
      dim = 1;
    else
      dim = 2;
    end
end

% handle empty S specially
if isempty(S)
  switch dim
    case 1
      meds = SegArray.constant(nan,1,sz(2));
    case 2
      meds = nan(sz(1),1);
    otherwise
      throwAsCaller(MException('Matlab:SegArray:DimTooHighMan',...
                  'Dimensions greater than 2 not supported by nanmedian.'));
  end
  return
end

% sort values (medians sort to end)
sorted = sort(S,dim);
mindex = (sum(~isnan(sorted),dim) + 1) / 2;
% calculate linear indices for median floor and ceiling
switch dim
  case 1
    col = 0:sz(1):numel(S)-1;
    uploc = col + ceil(mindex);
    dnloc = col + floor(mindex);
  case 2
    row = (1:sz(1))';
    uploc = (ceil(mindex)-1).*sz(1) + row;
    dnloc = (floor(mindex)-1).*sz(1) + row;
  otherwise
      throwAsCaller(MException('Matlab:SegArray:DimTooHighMan',...
                  'Dimensions greater than 2 not supported by nanmedian.'));
end
% look up median spanning values and average them
% (note: same location for odd numbers of non-nan values,
% two adjacent locations for even non-nans)
meds = ( subsref(sorted,substruct('()',{uploc})) + ...
         subsref(sorted,substruct('()',{dnloc})) ) / 2;

