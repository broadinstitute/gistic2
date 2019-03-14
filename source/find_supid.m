function idx=find_supid(D,supaccs,rc,varargin)
%FIND_SUPID Search the supacc field of data structure D for the identifier SUPACCS
%
%IDX = FIND_SUPID(D,SUPACCS,RC,VARARGIN)
%
%     RC: 'row' or 'col' (see is_col.m)
%
%     Revisions:
%           5 Oct 07:  Revise to allow for case mismatched supaccs.
%           22 Oct 07:  Put regexprep check to first argument of strmatch.
%           16 Dec 07:  Add support for supacc or gsupacc (default rc=row)
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


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
