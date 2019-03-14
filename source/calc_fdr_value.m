function fdr=calc_fdr_value(p)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if isempty(p)
    fdr = p;
    return
end

if size(p,1)==1
  p=p';
  trans=1;
else
  trans=0;
end

[sp,ord]=sort(p);

fdr=sp*size(p,1)./repmat((1:size(p,1))',1,size(p,2)); % make sure identical values
                                    % are taken care of 

fdr(fdr>1)=1;

fdr=[fdr; ones(1,size(fdr,2))];
for i=size(p,1):-1:1
  fdr(i,:)=min(fdr(i:(i+1),:),[],1);
end
fdr=fdr(1:(end-1),:);
ordmat=ord+repmat(0:size(p,1):size(p,1)*(size(fdr,2)-1),size(ord,1),1);
fdr(ordmat(:))=fdr(:);

if trans
  fdr=fdr';
end
