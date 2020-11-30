function cm = bluepink(N)
% bluepink - blue<=>pink colormap with true white at zero
%
%   CM = bluepink(N)
%
%   create a colormap with N levels, default N = 127
%
%

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

if ~exist('N','var')
    N = 127;
end

n = floor(N/2);
up = (0:(n-1))'/n;
dn = flipud(up);
hi = ones(n,1);
cm = [[up,up,hi];...
      [1, 1, 1];...
      [hi,dn,dn]];

if nargout==0
    colormap(cm);
    cm=[];
end
