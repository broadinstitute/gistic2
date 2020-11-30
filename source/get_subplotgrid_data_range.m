function r=get_subplotgrid_data_range(gr,i,j)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)


h=subplotgrid(gr,i,j,0);
% keyboard
r=caxis(h);


