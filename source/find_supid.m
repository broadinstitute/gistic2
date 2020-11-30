function idx = find_supid(D,supaccs,rc,varargin)
%FIND_SUPID Search the supacc field of data structure D for the identifier SUPACCS
%
% IDX = FIND_SUPID(D,SUPACCS,RC,VARARGIN)
%
%     RC: 'row' or 'col' (see is_col.m)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

if ischar(supaccs)
  supaccs=cellstr(supaccs);
end

if ~exist('rc','var')
  rc='cols';
end

idx=[];

if is_col(rc)
  if ~isempty(D.supacc)
    for i=1:length(supaccs)
      idx=[idx strmatch(lower(regexprep(supaccs{i},':.*','')),lower(regexprep(cellstr(D.supacc),':.*','')),varargin{:})];
    end
  else idx = [];
  end
else
  if ~isempty(D.gsupacc)
    for i=1:length(supaccs)
      idx=[idx strmatch(lower(regexprep(supaccs{i},':.*','')),lower(regexprep(cellstr(D.gsupacc),':.*','')),varargin{:})];
    end
  else idx = [];
  end
end
