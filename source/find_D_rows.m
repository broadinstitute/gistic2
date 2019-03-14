function idx=find_D_rows(D,gsupacc_name,val)
%FIND_D_ROWS Return indices of the rows in data structure that match the 
%specified field of .gsupacc with the given value.
%
%   IDX = FIND_D_ROWS(D,GSUPACC_NAME,VAL) returns indices IDX such that
%   D.dat(IDX,:) filters for rows whos GSUPACC_NAME matches VAL in the .gsupdat. 
%
%   Example: IDX = FIND_D_ROWS(D,'Interesting',1) returns the column
%   indices of the data rows with a value of 1 for 'Interesting' in .gsupdat 
%   field.  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


gsupid=strmatch(gsupacc_name,D.gsupacc,'exact');
if length(gsupid)~=1
  error([ mfilename ': non-unique or non-valid id']);
end
idx=find(D.gsupdat(gsupid,:)==val);
