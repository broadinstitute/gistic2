function [ext,fname_no_ext]=file_ext(fname)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

pos1=find(fname==filesep);
if isempty(pos1)
  st=1;
else
  st=pos1(end)+1;
end
if st>length(fname)
  error('fname is a directory');
end
  
pos=find(fname(st:end)=='.');
if ~isempty(pos)
  pos=pos+st-1;
  ext=fname((pos(end)+1):end);
  fname_no_ext=fname(1:(pos(end)-1));
else
  ext='';
  fname_no_ext=fname;
end
