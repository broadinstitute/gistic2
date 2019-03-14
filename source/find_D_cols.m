function idx=find_D_cols(D,supacc_name,val)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

supid=strmatch(supacc_name,regexprep(cellstr(D.supacc),':.*',''),'exact');
if length(supid)~=1
  error([ mfilename  ': non-unique or non-valid id']);
end
idx=find(D.supdat(supid,:)==val);
