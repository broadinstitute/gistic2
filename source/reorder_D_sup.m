function D=reorder_D_sup(D,rc,idx)
%REORDER_D_SUP reorder the supdat or the gsupdat's acc. 
%
%   D = reorder_D_sup(D,RC,IDX)
%           D = the data structure.
%           
%           RC = 'row' or 'col'.  Use 'row' to reorder gsupdat, 'col' to
%           reorder supdat.
%
%           IDX = reorder indices.

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

if isfield(D,'gsupdat') && strcmp(rc(1:3),'row')
  if isempty(idx)
    D=rmfield(D,{'gsupdat','gsupdesc','gsupacc'});
  else
    D.gsupdat=D.gsupdat(idx,:);
    D.gsupdesc=D.gsupdesc(idx,:);
    D.gsupacc=D.gsupacc(idx,:);
  end
  if isfield(D,'gsupmark')
    if isempty(idx)
      D=rmfield(D,{'gsupmark'});
    else
      D.supmark=D.gsupmark(idx);
    end
  end
end

if isfield(D,'supdat') && strcmp(rc(1:3),'col')
  if isempty(idx)
    D=rmfield(D,{'supdat','supdesc','supacc'});
  else
    D.supdat=D.supdat(idx,:);
    D.supdesc=D.supdesc(idx,:);
    D.supacc=D.supacc(idx,:);
  end
  if isfield(D,'supmark')
    if isempty(idx)
      D=rmfield(D,{'supmark'});
    else
      D.supmark=D.supmark(idx);
    end
  end
end

%this is a three-dimensional array with info from all merged arrays
if isfield(D,'mergsupdat') && strcmp(rc(1:3),'col')
  if isempty(idx)
    D=rmfield(D,{'mergsupdat','mergsupdesc','mergsupacc'});
  else
    D.mergsupdat=D.mergsupdat(idx,:,:);
    D.mergsupdesc=D.mergsupdesc(idx,:);
    D.mergsupacc=D.mergsupacc(idx,:);
  end
end



