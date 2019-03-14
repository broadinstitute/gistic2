function tab=read_table(fname,format,dlm,headerlines,varargin)
% read_table(fname,format,dlm,headerlines,varargin)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


fpos=0;
if headerlines~=0
  if ischar(fname)
    f=fopen(fname,'r');
  else
    f=fname;
    fpos=ftell(f);
  end
  if length(dlm)~=1
    tab.dlm=find_dlm(fname,dlm);
  else
    tab.dlm=dlm;
  end
  if headerlines>0
    tab.headers=read_dlm_file(f,dlm,headerlines);
  elseif headerlines==-1 % R convension
    headerlines=1;
    tab.headers=read_dlm_file(f,dlm,headerlines);
    tab.headers{1}=['EMPTY' tab.headers{1,:}];
  end
%   fclose(f);
else
  if ischar(fname)
    f=fopen(fname,'r');
  else
    f=fname;
    fpos=ftell(f);
  end
  tab.headers={};
  tab.dlm=dlm;
end

if isempty(format)
  if isempty(tab.headers)
    error('must have either format or headerlines');
  else
%    format=[repmat('%s',1,length(tab.headers{1})) '\n'];
    format=[repmat('%s',1,length(tab.headers{end})) '\n'];   % (to allow for multiple header lines)
  end
elseif iscell(format)
  if isempty(tab.headers)
    error('must have either format or headerlines');
  else
    if length(format)==1
%      format=[repmat(format{1},1,length(tab.headers{1})) '\n'];
      format=[repmat(format{1},1,length(tab.headers{end})) '\n'];
    else
%      [ length(tab.headers{1})-format{2} ]
%      format=[format{1} repmat(format{3},1,length(tab.headers{1})-format{2}) '\n'];
      format=[format{1} repmat(format{3},1,length(tab.headers{end})-format{2}) '\n'];
    end
  end
end
  

%format
%if strcmp(format(1:(end-2)),'\n')
if strcmp(format((end-1):end),'\n')      %%%%% fixed 2011-01-05 ML:  was causing incorrect handling of empty first column
  format=format(1:(end-2)); % remove '\n' at the end
end

verbose(['Reading file using format:' format],10);
fseek(f,fpos,'bof');
tab.dat=textscan(f,format,'headerLines',headerlines,'delimiter',tab.dlm,varargin{:});
fclose(f);

