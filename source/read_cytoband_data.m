function cyto=read_cytoband_data(fname)

%fname 'cytoBand_hg15.txt'

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

f=read_dlm_file(fname);


for i=1:length(f)
  cyto(i).chr=f{i}{1}(4:end);
  cyto(i).chrn=chromosome2num(cyto(i).chr);
  cyto(i).name=[ cyto(i).chr f{i}{4} ];
  cyto(i).start=str2num(f{i}{2});
  cyto(i).end=str2num(f{i}{3});
  cyto(i).stain=f{i}{5};
end
