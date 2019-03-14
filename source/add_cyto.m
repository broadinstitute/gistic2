function C=add_cyto(C,cyto)
% ADD_CYTO add cytoband information to copy number D-struct
%
%   D = add_cyto(D,CYTO)
%
% D is a copy number D-structure, CYTO is a struct array containing
% cytoband information. D is updated with new fields 'armn, 'cyton', and
% 'cyto_stain'.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if isfield(C,'cyton') && isfield(C,'armn')
  return;
end

if is_mb(C)
  error('not supporting Mb anymore');
end

if ~isfield(cyto,'chrn')
  x=chromosome2num({cyto(:).chr});
for i=1:length(x)
  cyto(i).chrn=x(i);
end
else
  for i=1:length({cyto.chrn})
    x(i)=cyto(i).chrn;
  end
end

C.cyton=zeros(length(C.pos),1);
C.armn=zeros(length(C.pos),1);
C.cyto_stain=nan(length(C.pos),1);
%C.cyto=cell(length(C.pos),1);
for i=1:max(x) % go over chromosomes
  fi=find(x==i); %find cytobands in each chromosome
  for j=1:length(fi)
    p=fi(j);
    st=cyto(p).start;
    en=cyto(p).end;
    snp=find(C.chrn==i & C.pos>=st & C.pos<en);
    if ~isempty(snp)
      C.cyton(snp)=p;
      %    C.cyto(snp)=cyto(p).name;
      if ~isempty(find(cyto(p).name=='p'))
        C.armn(snp)=1;
      else
        C.armn(snp)=2;
      end
      stain=cyto(p).stain;
      if ~strcmp(stain,'Empty')
        if length(stain)>4 && ~strcmp('stalk',stain)
          C.cyto_stain(snp)=str2num(stain(5:end));
        elseif length(stain)>=2 && stain(2)=='n'
          C.cyto_stain(snp)=0;
        end
      end
    end      
  end
  % fix for snps that are in the telomere
  snp=find(C.chrn==i & C.pos>cyto(fi(end)).end); % find snps beyond the last cytoband
  C.cyton(snp)=fi(end);
  C.armn(snp)=2;
end

