function write_array_list(fname,C)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

f=fopen(fname,'w');
fprintf(f,'%s\n','Array');
for i=1:size(C.dat,2)
  fprintf(f,'%s\r\n',deblank(char(C.sdesc(i,:))));
end
fclose(f);
