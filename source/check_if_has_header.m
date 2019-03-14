function [res,nfields]=check_if_has_header(fname,numeric_field_idx,varargin)
%[res,nfields] = check_if_has_header(fname,numeric_field_idx,varargin)
%res = 1 if has header
%nfields = number of columns in file

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


[fid,msg]=fopen(fname);
if fid==-1
    error ([msg,':',fname])
end
if exist('varargin','var') && ~isempty(varargin)
  ln=textscan(fgetl(fid),'%s',varargin{:});
else
  ln=textscan(fgetl(fid),'%s');
end
fclose(fid);
res=isempty(str2num(ln{1}{numeric_field_idx}));

nfields = length(ln);