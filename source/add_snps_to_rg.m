function rg=add_snps_to_rg(rg,C)
% adds snp locations to each gene
%
%---
% $Id$
% $Date: 2014-01-31 15:29:07 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%%% add chrn, snps, and nsnps for each gene to rg
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
