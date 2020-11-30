function u=convert(v,m)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

disp('FIXME');
[M,m1,m2]=match_string_sets_hash(num2str(m(:,1),'%f'),num2str(as_column(v),'%f'));

u=m(m1,2);
