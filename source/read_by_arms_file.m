function B=read_by_arms_file(fname)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)
    

  tbl = read_table(fname,{'%s',1,'%f'},char(9),1,'commentstyle','shell');
  B.sdesc = tbl.headers{1}(2:end);
  B.dat = cat(2,tbl.dat{2:end});
  B.armnames = tbl.dat{1};
