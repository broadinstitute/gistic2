function [res,resi]=grep(reg_exp,strs,res_is_idx)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)


if ischar(strs)
  strs=cellstr(strs);
end

resi=find(~cellfun('isempty',cat(1,regexp(strs,reg_exp))));
res=strs(resi);

if nargout==0
  disp(catn(res));
end

if exist('res_is_idx','var') && res_is_idx
  [res,resi]=exchange_vars(res,resi);
end
