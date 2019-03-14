function C=shrink_D(C)

% keep raw cbs_rl marker dat sdesc sis used_normals medians

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

C=rmfield_if_exists(C,{'origidx','supacc','supdat','supdesc','cbs_fixed','affy_calls','cbs','sm1','sm2','sm2j','sm3','smooth','temp','history'});
