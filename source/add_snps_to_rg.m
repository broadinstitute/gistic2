function rg = add_snps_to_rg(rg,C)
% adds snp locations to each gene

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

if ~isfield(rg,'chrn')
  rg = add_chrn(rg);
end

if ~isfield(rg,'snps')
  for i = 1:length(rg)
    rg(i).snps=find_snps(C,rg(i).chrn,rg(i).start,rg(i).end,2);
    rg(i).nsnps=length(rg(i).snps);
    rg(i).midsnp=floor(median(rg(i).snps));
  end
else
  verbose('Already has snps');
end
