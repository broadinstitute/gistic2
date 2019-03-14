function D=add_supmark(D,c)
% ADD_SUPMARK

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~exist('c','var')
  c=read_colorscheme();
end

for i=1:size(D.supdat,1)
  [st,sn]=break_sup_names(deblank(D.supacc(i,:)));
  m=max(ceil(max(D.supdat(i,:),[],2)),length(unique(D.supdat(i,:))));
  if ~isempty(sn) || m>2
    c1=repmat(c,max(length(sn),m),1);
    D.supmark(i).colormap=c1(1:max(length(sn),m),:);
  else
    D.supmark(i).colormap=[0 1 0; 1 0 0];
  end
  D.supmark(i).height=1;
  D.supmark(i).patchwidth=1;
  D.supmark(i).linewidth=0;
  D.supmark(i).marker=cellstr(repmat('.',max(length(sn),2),1));
  D.supmark(i).marker_size=5*ones(max(length(sn),2),1);
end
