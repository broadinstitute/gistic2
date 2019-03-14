function [st,en]=enlarge_range(C,st_snp,en_snp,before,after,min_st,max_en)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('after','var')
  after=before;
end

if ~isfield(C,'chrn')
  C=add_chrn(C);
end

chr=C.chrn(st_snp);
chr_en=C.chrn(en_snp);

if chr~=chr_en
  error('Start and End are not on the same chromosome');
end

in_chr=find(C.chrn==chr);
begin_chr=min(in_chr);
end_chr=max(in_chr);

st=max(begin_chr,st_snp-before);
en=min(end_chr,en_snp+after);

if exist('min_st','var')
  st=max(min_st,st);
end

if exist('max_en','var')
  en=min(max_en,en);
end


