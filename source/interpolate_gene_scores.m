function score = interpolate_gene_scores(D,gg_rg,field,interp_fun)
%INTERPOLATE_GENE_SCORES marker representation for gene gistic scores
%
%   SCORE = interpolate_gene_scores(D,GG_RG,FIELD,INTERP_FUN)
%
% Returns a by-marker projection of gene values in GG_RG. Values between
% two genes on a chromosome are interpolated. The D structure provides
% the mapping of the genome to the markers. FIELD names the field in
% GG_RG to project (e.g. 'ads' or 'q') and INTERP_FUN is a function
% handle to the interpolating function (e.g. @max or @min).

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% project gene scores onto blank (NaN) genomic markers
q2 = nan(size(D.pos));
for g = 1:length(gg_rg)
  gsnps = gg_rg(g).snps;
  q2(gsnps) = interp_fun(q2(gsnps),gg_rg(g).(field));
end

% segment the data into a runlength and add a 4th column of chromosomes
rlq2 = runlength(q2,runlength(D.chrn));
rlq2 = [rlq2,D.chrn(rlq2(:,1))];

% NaNs at ends of chromosomes are given value of next non-NaN inward
chrends = [find(diff(rlq2(:,4)));size(rlq2,1)];
nanatend = chrends(isnan(rlq2(chrends,3)));
rlq2(nanatend,3) = rlq2(nanatend-1,3);
chrbegs = [1;chrends(1:end-1)+1];
nanatbeg = chrbegs(isnan(rlq2(chrbegs,3)));
rlq2(nanatbeg,3) = rlq2(nanatbeg+1,3);
% remaining NaN markers are interpolated
naninmid = find(isnan(rlq2(:,3)));
rlq2(naninmid,3) = interp_fun(rlq2(naninmid+1,3),rlq2(naninmid-1,3));
% rebuild the full markers
score = derunlength(rlq2);
