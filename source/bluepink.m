function cm = bluepink(N)
% bluepink - blue<=>pink colormap with true white at zero
%
%   CM = bluepink(N)
%
%   create a colormap with N levels, default N = 127
%
%

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
