function B = read_by_gene_file(fname,refgene)
%READ_BY_GENE_FILE

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

tbl = read_table(fname,{'%s%s%s',3,'%f'},char(9),1,'commentstyle','shell');
B.sdesc = tbl.headers{1}(4:end);
B.dat = cat(2,tbl.dat{4:end});
B.gacc = tbl.dat{1};
B.gdesc = tbl.dat{2};
B.cyto = tbl.dat{3};

if exist('refgene','var')
  refgene = load_refgene(refgene);
  rg = refgene.gene;
  
  [~,m1,m2] = match_string_sets_hash({rg.symb},B.gacc);
  [~,umi,~] = lunique(m2);
  um1 = m1(umi);
  B.chr = {rg(um1).chr};
  B.chrn = chromosome2num(B.chr);
  B.pos = floor(0.5*(cat(1,rg(um1).start)+cat(1,rg(um1).end)));
end
