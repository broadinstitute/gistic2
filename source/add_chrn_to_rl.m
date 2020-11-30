function C=add_chrn_to_rl(C)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)


for i=1:length(C.cbs_rl)
  cur=C.cbs_rl{i};
  cur=[cur zeros(size(cur,1),1)];
  for j=1:size(cur,1)
%    disp([ i j cur(j,1:2)]);
    cur(j,end)=C.chrn(cur(j,1));
  end
  C.cbs_rl{i}=cur;
end
