function D=load_D(fname,varname)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if exist('varname','var')
  D=load(fname,varname);
else
  tmp=load(fname);
end

nms=fieldnames(tmp);
D=getfield(tmp,nms{1});

if isfield(D,'marker') && ischar(D.marker)
  D.marker=cellstr(D.marker);
end

if isfield(D,'chr') && ischar(D.chr)
  D.chr=cellstr(D.chr);
end

if isfield(D,'sdesc') && ischar(D.sdesc)
  D.sdesc=cellstr(D.sdesc);
end
